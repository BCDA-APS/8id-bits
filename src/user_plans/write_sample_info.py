import yaml
from pathlib import Path


USER_PLAN_DIR = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans")
SAMPLE_INFO_FILE = USER_PLAN_DIR / "sample_info.yaml"


disp = 0.7  # Granite is at 923.0 as of 03/18/2026

y_center = 28.66

x_radius = 0.5
y_radius = 0.5

x_points = 25
y_points = 25

inner_range = 2 * x_radius
outer_range = 2 * y_radius

x_centers = [
    303.0,
    298.0,
    293.2,
    269.7,
    264.8,
    259.9,
    236.5,
    231.5,
    226.5,
    188.7,
    183.7,
    178.7,
    155.4,
    150.4,
    145.4,
    122.1,
    117.1,
    112.1,
    74.6,
    69.6,
    64.6,
    41.4,
    36.4,
    31.4,
    8.2,
    3.2,
]

headers = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


sample_info = {
    "defaults": {
        "inner_motor": "sample.x",
        "outer_motor": "sample.y",
        "inner_range": inner_range,
        "outer_range": outer_range,
        "inner_pts": x_points,
        "outer_pts": y_points,
    },
    "samples": {},
}


for i in range(1, 27):
    sample_key = f"sample_{i}"

    sample_info["samples"][sample_key] = {
        "sample_name": f"Test{i}",
        "header": headers[i - 1],
        "inner_center": x_centers[i - 1] + disp,
        "outer_center": y_center,
    }


def write_sample_info():
    USER_PLAN_DIR.mkdir(parents=True, exist_ok=True)

    yaml_text = yaml.safe_dump(sample_info, sort_keys=False)
    yaml_text = yaml_text.replace("\n  sample_", "\n\n  sample_")

    with open(SAMPLE_INFO_FILE, "w") as f:
        f.write(yaml_text)

    print(f"Wrote {SAMPLE_INFO_FILE}")


if __name__ == "__main__":
    write_sample_info()