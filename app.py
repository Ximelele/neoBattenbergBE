import os

from flask import Flask, jsonify, request
from flask_cors import CORS

from Routes.scariance_api import scariance


app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "http://localhost:5173"}})

app.config['IMAGE_DIRECTORY'] = os.environ.get('IMAGE_DIRECTORY', '/default/images')
app.config['DATA_DIRECTORY'] = os.environ.get('DATA_DIRECTORY', '/default/data')
app.register_blueprint(scariance)

if __name__ == '__main__':
    app.run(debug=True, port=6800)
