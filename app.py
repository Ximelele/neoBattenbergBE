from flask import Flask, jsonify, request
from flask_cors import CORS  # Import flask-cors

from Anotations.StrVCTVRE_annot import start_analysis
from Routes.annotSV_api import annotsv_BP
from Routes.battenberg_api import battenberg_BP
from Routes.strvctre_api import strvctrvreBP

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "http://localhost:5173"}})
app.register_blueprint(strvctrvreBP)
app.register_blueprint(annotsv_BP)
app.register_blueprint(battenberg_BP)


@app.route('/api/update-threshold', methods=['POST'])
def bar_plot():
    # Get JSON data from the request
    data = request.get_json()
    if not data or 'threshold' not in data:
        return jsonify({"error": "Invalid input, 'threshold' is required"}), 400
    print(data['threshold'])
    # Call your function and return the result
    result = start_analysis("/Users/martindruzbacky/PycharmProjects/DP_helper/patient1_annotated.vcf",
                            data['threshold'])
    return jsonify(result)


if __name__ == '__main__':
    app.run(debug=True, port=6800)
