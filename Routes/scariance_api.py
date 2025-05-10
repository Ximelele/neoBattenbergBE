import pandas as pd
from flask import jsonify, Blueprint, send_from_directory, request,current_app
import os
import re

scariance = Blueprint('scariance', __name__)


@scariance.route('/api/images', methods=['GET'])
def list_images():
    if not os.path.exists(current_app.config['IMAGE_DIRECTORY']):
        return jsonify({"error": "Image directory does not exist"}), 404

    try:
        # Get all image files in the directory
        image_files = [f for f in os.listdir(current_app.config['IMAGE_DIRECTORY']) if f.endswith(('.png', '.jpg', '.jpeg', '.svg'))]

        def extract_chr(file_name):
            match = re.search(r'chr(\d+)', file_name)
            return int(match.group(1)) if match else float('inf')

        image_files.sort(key=extract_chr)

        categorized_images = {
            "heterozygousData": [],
            "RAFseg": [],
            "segment": [],
            "copynumberprofile": [],
            "cnv_per_chrom": [],
        }

        for file in image_files:
            url = f"http://127.0.0.1:6800/{file}"

            if "heterozygousData" in file:
                categorized_images["heterozygousData"].append(url)
            elif "RAFseg" in file:
                categorized_images["RAFseg"].append(url)
            elif "segment" in file:
                categorized_images["segment"].append(url)
            elif any(key in file for key in ["copynumberprofile", "distance", "nonroundedprofile"]):
                categorized_images["copynumberprofile"].append(url)
            elif "cnv_per_chrom" in file:
                categorized_images["cnv_per_chrom"].append(url)

        return jsonify(categorized_images)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@scariance.route("/api/copy-number-data", methods=["GET"])
def get_copy_number_data():
    """Return the full copy‑number table as JSON (no pagination)."""
    try:
        # Absolute path to the copy‑number file
        file_path = os.path.join(
            current_app.config['DATA_DIRECTORY']
        )

        # Fail fast if the file is missing
        if not os.path.exists(file_path):
            return jsonify({"error": "Copy number data file not found"}), 404

        # Load the *entire* tab‑delimited file
        df = pd.read_csv(file_path, sep="\t")

        # Build the response without any paging metadata
        response = {
            "columns": df.columns.tolist(),
            "data": df.to_dict(orient="records"),
        }

        return jsonify(response)

    except Exception as e:
        # Catch‑all for unexpected errors
        return jsonify({"error": str(e)}), 500



@scariance.route('/<path:filename>', methods=['GET'])
def serve_image(filename):
    """
    Serve individual image files from the directory.
    """
    return send_from_directory(current_app.config['IMAGE_DIRECTORY'], filename)
