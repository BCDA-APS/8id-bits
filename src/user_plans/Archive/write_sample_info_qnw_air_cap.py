import yaml

disp = 0.44  # Granite is at 923.0 as of 03/18/2026
y_cen = 20.84
x_radius = 0.5
y_radius = 0.5
x_pts = 25 
y_pts = 25

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
        "header": "R",
        "temp_zone": "qnw_env1"
    },
    "sample_1": {
        "sample_name": "Test1",
        "x_cen": 303.0+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Aa",
        "temp_zone": "qnw_env1"
    }, 
    "sample_2": {
        "sample_name": "Test2",
        "x_cen": 298.0+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ab",
        "temp_zone": "qnw_env1"
    },
    "sample_3": {
        "sample_name": "Test3",
        "x_cen": 293.2+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ac",
        "temp_zone": "qnw_env1"
    },
    "sample_4": {
        "sample_name": "Test4",
        "x_cen": 269.7+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ba",
        "temp_zone": "qnw_env1"
    },
    "sample_5": {
        "sample_name": "Test5",
        "x_cen": 264.8+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Bb",
        "temp_zone": "qnw_env1"
    },
    "sample_6": {
        "sample_name": "Test6",
        "x_cen": 259.9+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Bc",
        "temp_zone": "qnw_env1"
    },
    "sample_7": {
        "sample_name": "Test7",
        "x_cen": 236.5+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ca",
        "temp_zone": "qnw_env1"
    },
    "sample_8": {
        "sample_name": "Test8",
        "x_cen": 231.5+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Cb",
        "temp_zone": "qnw_env1"
    },
    "sample_9": {
        "sample_name": "Test9",
        "x_cen": 226.5+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Cc",
        "temp_zone": "qnw_env1"
    },
    "sample_10": {
        "sample_name": "Test10",
        "x_cen": 188.7+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Da",
        "temp_zone": "qnw_env2"
    },
    "sample_11": {
        "sample_name": "Test11",
        "x_cen": 183.7+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Db",
        "temp_zone": "qnw_env2"
    },
    "sample_12": {
        "sample_name": "Test12",
        "x_cen": 178.7+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Dc",
        "temp_zone": "qnw_env2"
    },
    "sample_13": {
        "sample_name": "Test13",
        "x_cen": 155.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ea",
        "temp_zone": "qnw_env2"
    },
    "sample_14": {
        "sample_name": "Test14",
        "x_cen": 150.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Eb",
        "temp_zone": "qnw_env2"
    },
    "sample_15": {
        "sample_name": "Test15",
        "x_cen": 145.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ec",
        "temp_zone": "qnw_env2"
    },
    "sample_16": {
        "sample_name": "Test16",
        "x_cen": 122.1+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Fa",
        "temp_zone": "qnw_env2"
    },
    "sample_17": {
        "sample_name": "Test17",
        "x_cen": 117.1+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Fb",
        "temp_zone": "qnw_env2"
    },
    "sample_18": {
        "sample_name": "Test18",
        "x_cen": 112.1+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Fc",
        "temp_zone": "qnw_env2"
    },
    "sample_19": {
        "sample_name": "Test19",
        "x_cen": 74.6+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ga",
        "temp_zone": "qnw_env3"
    },
    "sample_20": {
        "sample_name": "Test20",
        "x_cen": 69.6+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Gb",
        "temp_zone": "qnw_env3"
    },
    "sample_21": {
        "sample_name": "Test21",
        "x_cen": 64.6+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Gc",
        "temp_zone": "qnw_env3"
    },
    "sample_22": {
        "sample_name": "Test22",
        "x_cen": 41.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ha",
        "temp_zone": "qnw_env3"
    },
    "sample_23": {
        "sample_name": "Test23",
        "x_cen": 36.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Hb",
        "temp_zone": "qnw_env3"
    },
    "sample_24": {
        "sample_name": "Test24",
        "x_cen": 31.4+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Hc",
        "temp_zone": "qnw_env3"
    },
    "sample_25": {
        "sample_name": "Test25",
        "x_cen": 8.2+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ia",
        "temp_zone": "qnw_env3"
    },
    "sample_26": {
        "sample_name": "Test26",
        "x_cen": 3.2+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ib",
        "temp_zone": "qnw_env3"
    },
    "sample_27": {
        "sample_name": "Test27",
        "x_cen": -1.9+disp,
        "y_cen": y_cen,
        "x_radius": x_radius,
        "y_radius": y_radius,
        "x_pts": x_pts,
        "y_pts": y_pts,
        "header": "Ic",
        "temp_zone": "qnw_env3"
    }
}

# Writing to a file
# Generate the YAML string first
yaml_str = yaml.dump(nested_dict, default_flow_style=False, sort_keys=False)

# Add empty lines between samples.
# We replace a newline followed by 'sample_' with two newlines followed by 'sample_'.
# This works safely because top-level keys ('sample_0', etc.) are not indented, 
# while nested keys (e.g., 'sample_name') would be indented and won't match.
formatted_yaml = yaml_str.replace('\nsample_', '\n\nsample_')

with open("sample_info.yaml", "w") as f:
    f.write(formatted_yaml)