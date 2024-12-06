from flask import jsonify, Blueprint

from Models.FileLoad import load_patient_data

strvctrvreBP = Blueprint('strvctre', __name__)


@strvctrvreBP.route('/api/patient', methods=['GET'])
def get_strvctrvre_patient_data():
    try:
        patient_data = load_patient_data('/Users/martindruzbacky/PycharmProjects/DP_helper/strvctrvre.json')

        return jsonify({"patient_data": patient_data}), 200
    except FileNotFoundError as e:
        return jsonify({"error": str(e)}), 404
    except Exception as e:
        return jsonify({"error": str(e)}), 500
