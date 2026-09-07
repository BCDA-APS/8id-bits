import json

disp = 10.0  # Granite is at 923.0 as of 04/22/2025
y_cen = 21.15
# Sample nested dictionary
nested_dict = {
    "sample_0": {
        "sample_name": "TestRheo",
        "x_cen": -0.65,
        "y_cen": 349.783,
        "x_radius": 0.0,
        "y_radius": 0.0,
        "x_pts": 20,
        "y_pts": 20,
        "header": "X",
        "temp_zone": "qnw_env1"
    },
    "sample_1": {
        "sample_name": "Test1",
        "x_cen": 287.3 + disp, 
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "A",
        "temp_zone": "qnw_env1"
    },
    "sample_2": {
        "sample_name": "Test2",
        "x_cen": 254.4 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "B",
        "temp_zone": "qnw_env1"
    },
    "sample_3": {
        "sample_name": "Test3",
        "x_cen": 221.3 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "C",
        "temp_zone": "qnw_env1"
    },
    "sample_4": {
        "sample_name": "Test4",
        "x_cen": 173.3 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "D",
        "temp_zone": "qnw_env2"
    },
    "sample_5": {
        "sample_name": "Test5",
        "x_cen": 140.0 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "E",
        "temp_zone": "qnw_env2"
    },
    "sample_6": {
        "sample_name": "Test6",
        "x_cen": 107.3 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "F",
        "temp_zone": "qnw_env2"
    },
    "sample_7": {
        "sample_name": "Test7",
        "x_cen": 59.3 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "G",
        "temp_zone": "qnw_env3"
    },
    "sample_8": {
        "sample_name": "Test8",
        "x_cen": 26.3 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "H",
        "temp_zone": "qnw_env3"
    },
    "sample_9": {
        "sample_name": "Test9",
        "x_cen": -7.7 + disp,
        "y_cen": y_cen,
        "x_radius": 0.5,
        "y_radius": 0.5,
        "x_pts": 20,
        "y_pts": 20,
        "header": "I",
        "temp_zone": "qnw_env3"
    }
}

# Writing to a file
with open("sample_info.json", "w") as f:
    json.dump(nested_dict, f, indent=4)


