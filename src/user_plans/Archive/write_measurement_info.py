import json
import re
import copy 

def write_measurement_info(sam_range):

    file_path = "/home/beams10/8IDIUSER/bluesky/src/user_plans/"

    all_samples_data = {}

    file_template = 'measurement_info_framework.json'
    with open(file_path+file_template, 'r') as f:
        base_sample_structure = json.load(f)

    # Loop to create sample_0 to sample_27
    for i in sam_range:
        sample_key = f"sample_{i}"
        all_samples_data[sample_key] = copy.deepcopy(base_sample_structure) # Use deepcopy to avoid modifying the same object

    # Step 1: Recursively convert all lists in the data to compact JSON strings.
    def convert_lists_to_compact_strings(obj):
        if isinstance(obj, dict):
            return {k: convert_lists_to_compact_strings(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            # Convert the list to its compact JSON string representation.
            # separators=(',', ':') removes any whitespace after commas and colons.
            return json.dumps(obj, separators=(',', ':'))
        else:
            return obj

    # Apply the conversion to the entire dataset
    processed_data = convert_lists_to_compact_strings(all_samples_data)

    # Step 2: Dump the processed data to a string with indentation.
    # This will put keys on new lines, but our list-strings will be quoted.
    # We use indent=2 to match the example's common indentation.
    json_output_string = json.dumps(processed_data, indent=2)

    # Step 3: Post-process the JSON string to remove the quotes around
    # the previously stringified lists.
    # This regex replaces a pattern like "[[...]]" with just [[...]]
    final_json_string = re.sub(r'"(\[.*?\])"', r'\1', json_output_string)

    file_name="measurement_info.json"
    with open(file_path+file_name, 'w') as f:
        f.write(final_json_string)
    
    print(f"Successfully generated '{file_name}' with {len(sam_range)} samples and custom formatting.")

        
