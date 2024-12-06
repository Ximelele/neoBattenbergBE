import os
import json


def load_patient_data(path: str):
    # Ensure the file exists
    file_path = path
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    with open(file_path, 'r') as file:
        data = file.read()
        # Check if the file content is JSON
        try:
            return json.loads(data)  # Parse JSON content
        except json.JSONDecodeError:
            return data
