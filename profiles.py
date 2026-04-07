"""
profiles.py

Load and apply factor schemas and vehicle profiles for terrain mobility
normalization.

The normalization chain:
    Raw cell value (from Surface Model GDB)
        → Factor Schema (LANDFIRE, USCS, etc.): code → terrain class
        → Vehicle Profile: terrain class → normalized mobility (0.0–1.0)

For continuous factors (F1 Elevation, F2 Slope), the factor schema is
skipped and the vehicle profile directly normalizes the raw value using
min/max/threshold parameters.

For categorical factors (F3 Vegetation, F4 Soil), the chain is:
    raw integer/string → schema lookup → "class" field → vehicle terrain_mobility dict

Usage:
    from profiles import ProfileEngine

    engine = ProfileEngine(
        factor_schemas=["profiles/factors/landfire_gap_2011.json",
                        "profiles/factors/uscs_soils.json"],
        vehicle_path="profiles/vehicles/light_wheeled.json",
    )

    # Normalize a single cell's raw values
    raw = {"F1": 2500.0, "F2": 12.0, "F3": 151, "F4": "CL"}
    normalized = engine.normalize(raw)
    # -> {"F1": 0.44, "F2": 0.66, "F3": 0.25, "F4": 0.55}

    # Compute mobility from normalized values
    mobility = engine.mobility(normalized)
    # -> (F1 / F2) * F3 * F4
"""

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


# Data Structures

@dataclass
class FactorSchema:
    """A loaded factor schema (e.g., LANDFIRE vegetation, USCS soils)."""

    schema_id: str
    factor: int               # 1-4
    factor_name: str
    code_type: str            # "integer" or "string"
    lookup: Dict[str, Dict[str, Any]]   # code -> {"class": ..., "label": ...}
    default: Dict[str, Any]   # Fallback for unknown codes

    def resolve(self, raw_value: Any) -> str:
        """
        Map a raw raster value to a terrain class string.

        Handles type coercion: integer codes are looked up as strings,
        NaN/None returns the default class.
        """
        if raw_value is None:
            return self.default.get("class", "unknown")

        # Handle NaN
        if isinstance(raw_value, float) and math.isnan(raw_value):
            return self.default.get("class", "unknown")

        # Coerce to lookup key
        key = str(int(raw_value)) if self.code_type == "integer" else str(raw_value).strip()

        entry = self.lookup.get(key)
        if entry is None:
            return self.default.get("class", "unknown")

        return entry.get("class", "unknown")


@dataclass
class ContinuousNorm:
    """Normalization parameters for a continuous factor (F1, F2)."""

    norm_min: float = 0.0
    norm_max: float = 100.0
    impassable_above: Optional[float] = None
    invert: bool = True  # True = lower raw value -> higher mobility

    def normalize(self, raw_value: float) -> float:
        """
        Normalize a raw continuous value to 0.0–1.0.

        If invert=True (default), lower raw values produce higher mobility:
            norm = 1.0 - (val - min) / (max - min)

        Values above impassable_above return 0.0.
        """
        if raw_value is None or (isinstance(raw_value, float) and math.isnan(raw_value)):
            return 0.0

        val = float(raw_value)

        # Impassable gate
        if self.impassable_above is not None and val > self.impassable_above:
            return 0.0

        # Clamp to range
        val = max(self.norm_min, min(val, self.norm_max))

        # Normalize
        span = self.norm_max - self.norm_min
        if span <= 0:
            return 1.0

        if self.invert:
            return 1.0 - (val - self.norm_min) / span
        else:
            return (val - self.norm_min) / span


@dataclass
class CategoricalNorm:
    """Normalization parameters for a categorical factor (F3, F4)."""

    terrain_mobility: Dict[str, float]   # terrain_class -> mobility 0.0-1.0
    default_mobility: float = 0.4

    def normalize(self, terrain_class: str) -> float:
        """Look up the mobility value for a terrain class."""
        return self.terrain_mobility.get(terrain_class, self.default_mobility)


@dataclass
class VehicleProfile:
    """A loaded vehicle profile with factor-specific normalization."""

    vehicle_id: str
    vehicle_name: str
    max_speed_mph: float
    max_speed_ms: float = 0.0  # Computed from mph
    description: str = ""

    # Per-factor normalization: factor key (F1-F4) -> normalizer
    continuous_norms: Dict[str, ContinuousNorm] = field(default_factory=dict)
    categorical_norms: Dict[str, CategoricalNorm] = field(default_factory=dict)

    def __post_init__(self):
        self.max_speed_ms = self.max_speed_mph * 0.44704


# Loading Functions

def load_factor_schema(path: Union[str, Path]) -> FactorSchema:
    """Load a factor schema JSON file."""
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    return FactorSchema(
        schema_id=data.get("schema_id", path.stem),
        factor=int(data["factor"]),
        factor_name=data.get("factor_name", f"Factor {data['factor']}"),
        code_type=data.get("code_type", "integer"),
        lookup=data.get("lookup", {}),
        default=data.get("default", {"class": "unknown"}),
    )


def load_vehicle_profile(path: Union[str, Path]) -> VehicleProfile:
    """Load a vehicle profile JSON file."""
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    profile = VehicleProfile(
        vehicle_id=data.get("vehicle_id", path.stem),
        vehicle_name=data.get("vehicle_name", path.stem),
        max_speed_mph=float(data.get("max_speed_mph", 45)),
        description=data.get("description", ""),
    )

    factors = data.get("factors", {})

    for factor_key, spec in factors.items():
        factor_type = spec.get("type", "continuous")

        if factor_type == "continuous":
            profile.continuous_norms[factor_key] = ContinuousNorm(
                norm_min=float(spec.get("min", 0)),
                norm_max=float(spec.get("max", 100)),
                impassable_above=spec.get("impassable_above"),
                invert=spec.get("normalization", "inverse_linear") == "inverse_linear",
            )

        elif factor_type == "categorical":
            terrain_mob = spec.get("terrain_mobility", {})
            # Convert all values to float
            terrain_mob = {k: float(v) for k, v in terrain_mob.items()}
            profile.categorical_norms[factor_key] = CategoricalNorm(
                terrain_mobility=terrain_mob,
                default_mobility=terrain_mob.get("unknown", 0.4),
            )

    return profile


# Profile Engine

# Minimum mobility threshold — below this, a cell is impassable
IMPASSABLE_THRESHOLD = 0.01


class ProfileEngine:
    """
    Chains factor schemas and vehicle profiles to normalize raw cell values
    into mobility scores.

    The engine handles two paths:
      - Continuous factors (F1, F2): raw value → vehicle continuous norm → 0.0-1.0
      - Categorical factors (F3, F4): raw value → schema lookup → terrain class
                                       → vehicle categorical norm → 0.0-1.0
    """

    def __init__(
        self,
        factor_schemas: Optional[List[Union[str, Path]]] = None,
        vehicle_path: Optional[Union[str, Path]] = None,
    ):
        # Factor schemas keyed by factor number (1-4)
        self.schemas: Dict[int, FactorSchema] = {}

        if factor_schemas:
            for schema_path in factor_schemas:
                schema = load_factor_schema(schema_path)
                self.schemas[schema.factor] = schema
                print(f"  Loaded factor schema: {schema.schema_id} "
                      f"(F{schema.factor} {schema.factor_name}, "
                      f"{len(schema.lookup)} codes)")

        # Vehicle profile
        if vehicle_path:
            self.vehicle = load_vehicle_profile(vehicle_path)
            print(f"  Loaded vehicle profile: {self.vehicle.vehicle_name} "
                  f"({self.vehicle.max_speed_mph} mph)")
        else:
            # Default generic vehicle
            self.vehicle = VehicleProfile(
                vehicle_id="generic",
                vehicle_name="Generic",
                max_speed_mph=30,
                description="Generic vehicle with no schema-based normalization",
            )

    def normalize_cell(self, raw_values: Dict[str, Any]) -> Dict[str, float]:
        """
        Normalize all factor values for a single cell.

        Parameters
        ----------
        raw_values : dict
            Raw values keyed by factor column name (F1, F2, F3, F4).
            Values can be float (continuous), int (categorical code),
            str (categorical code), or None/NaN.

        Returns
        -------
        dict[str, float]
            Normalized values keyed by factor column name, each in [0.0, 1.0].
        """
        normalized = {}

        for factor_key in ["F1", "F2", "F3", "F4"]:
            raw = raw_values.get(factor_key)

            # Check if this factor has a continuous normalizer
            if factor_key in self.vehicle.continuous_norms:
                norm = self.vehicle.continuous_norms[factor_key]
                normalized[factor_key] = norm.normalize(raw)

            # Check if this factor has a categorical normalizer
            elif factor_key in self.vehicle.categorical_norms:
                cat_norm = self.vehicle.categorical_norms[factor_key]

                # Get the factor number to find the right schema
                factor_num = int(factor_key[1])
                schema = self.schemas.get(factor_num)

                if schema is not None:
                    # Chain: raw -> schema -> terrain class -> vehicle mobility
                    terrain_class = schema.resolve(raw)
                    normalized[factor_key] = cat_norm.normalize(terrain_class)
                else:
                    # No schema — try direct lookup (raw value IS the class)
                    if raw is not None and not (isinstance(raw, float) and math.isnan(raw)):
                        terrain_class = str(raw).strip()
                        normalized[factor_key] = cat_norm.normalize(terrain_class)
                    else:
                        normalized[factor_key] = 0.0

            else:
                # No normalizer defined for this factor — default to 1.0 (no restriction)
                normalized[factor_key] = 1.0 if raw is None else 1.0

        return normalized

    def mobility(self, normalized: Dict[str, float]) -> float:
        """
        Compute mobility from normalized factor values.

        Mobility = (F1 / F2) * F3 * F4

        All inputs are 0.0–1.0. F2 is floored at IMPASSABLE_THRESHOLD
        to avoid division by zero.

        Returns a value where higher = faster traversal. Can exceed 1.0
        when F1 > F2.
        """
        f1 = normalized.get("F1", 1.0)
        f2 = normalized.get("F2", 1.0)
        f3 = normalized.get("F3", 1.0)
        f4 = normalized.get("F4", 1.0)

        # Any factor at zero = impassable
        if (f1 <= IMPASSABLE_THRESHOLD or f3 <= IMPASSABLE_THRESHOLD
                or f4 <= IMPASSABLE_THRESHOLD):
            return 0.0

        f2_safe = max(f2, IMPASSABLE_THRESHOLD)

        return max(0.0, (f1 / f2_safe) * f3 * f4)

    def cell_speed_ms(self, raw_values: Dict[str, Any]) -> float:
        """
        Compute effective traversal speed (m/s) for a cell from raw values.

        Full chain: raw values → normalize → mobility → speed.
        """
        normalized = self.normalize_cell(raw_values)
        mob = self.mobility(normalized)

        if mob <= IMPASSABLE_THRESHOLD:
            return 0.0

        # Cap mobility at 1.0 for speed calculation
        return self.vehicle.max_speed_ms * min(mob, 1.0)

    def is_impassable(self, raw_values: Dict[str, Any]) -> bool:
        """Check if a cell is impassable given its raw factor values."""
        normalized = self.normalize_cell(raw_values)
        return self.mobility(normalized) <= IMPASSABLE_THRESHOLD