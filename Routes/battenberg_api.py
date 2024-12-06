from flask import jsonify, Blueprint, send_from_directory
import os
import re

IMAGE_DIRECTORY = "/Users/martindruzbacky/PycharmProjects/DP_helper/test"

battenberg_BP = Blueprint('battenberg_BP', __name__)


@battenberg_BP.route('/api/images', methods=['GET'])
def list_images():
    """
    List all images in the specified directory and return their URLs.
    """
    if not os.path.exists(IMAGE_DIRECTORY):
        return jsonify({"error": "Image directory does not exist"}), 404

    try:
        # Get all image files in the directory
        image_files = [f for f in os.listdir(IMAGE_DIRECTORY) if f.endswith(('.png', '.jpg', '.jpeg', '.svg'))]

        def extract_chr(file_name):
            match = re.search(r'chr(\d+)', file_name)  # Look for 'chr<number>' in the file name
            return int(match.group(1)) if match else float('inf')  # Use infinity for files without 'chr'

        # Sort files by chromosome number
        image_files.sort(key=extract_chr)

        # Generate URLs for each image
        image_urls = [f"http://127.0.0.1:6800/{file}" for file in image_files]
        return jsonify({"images": image_urls})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@battenberg_BP.route('/<path:filename>', methods=['GET'])
def serve_image(filename):
    """
    Serve individual image files from the directory.
    """
    return send_from_directory(IMAGE_DIRECTORY, filename)
