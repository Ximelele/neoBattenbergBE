from Models.FileLoad import load_patient_data
from flask import  jsonify,Blueprint
annotsv_BP = Blueprint('annotsv_BP', __name__)

@annotsv_BP.route('/api/patient-annotsv', methods=['GET'])
def get_annotsv_patient_data():
    try:
        patient_data = load_patient_data('/Users/martindruzbacky/PycharmProjects/DP_helper/annotsv.json')
        return jsonify({"patient_data": patient_data}), 200
    except FileNotFoundError:
        return jsonify({"error": "Patient data file not found.", "code": 404}), 404
    except Exception as e:
        return jsonify({"error": f"Unexpected server error: {str(e)}", "code": 500}), 500