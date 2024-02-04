from flask import Flask, request, jsonify
import sqlite3
from server.refine import *
app = Flask(__name__)


@app.after_request
def add_cors_headers(response):
    response.headers.add('Access-Control-Allow-Origin', 'http://192.168.60.116:3000')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type, x-requested-with')
    return response


@app.route('/', methods=['GET'])
def hello_world():
    return 'Hello, World!'


@app.route('/api/add', methods=['POST'])
def add_numbers():
    data = request.get_json()
    num1 = data['num1']
    num2 = data['num2']
    result = num1 + num2
    return jsonify({'result': result})


@app.route('/datasets', methods=['GET'])
def init_dataset():
    """

    """

    # 连接到数据库
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # 执行PRAGMA命令来获取表的属性名
    cursor.execute("PRAGMA table_info(datasets)")

    # 获取结果并输出属性名
    attributes = [row[1] for row in cursor.fetchall()]

    # 执行查询语句，获取所有数据
    cursor.execute("SELECT * FROM datasets")

    # 获取查询结果
    rows = cursor.fetchall()

    response = jsonify({'attributes': attributes, 'data': rows})

    # 关闭游标和连接
    cursor.close()
    conn.close()

    return response


@app.route('/strategies', methods=['GET'])
def get_strategies():
    # 连接到数据库
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # 执行PRAGMA命令来获取表的属性名
    cursor.execute("PRAGMA table_info(technology_table)")

    # 获取结果并输出属性名
    attributes = [row[1] for row in cursor.fetchall()]

    # 执行查询语句，获取所有数据
    cursor.execute("SELECT * FROM technology_table")

    # 获取查询结果
    rows = cursor.fetchall()

    response = jsonify({'attributes': attributes, 'data': rows})

    # 关闭游标和连接
    cursor.close()
    conn.close()

    return response


@app.route('/refine', methods=['POST'])
def refine():
    Refiners = ['Eager Refiner', 'Lazy Refiner', 'LabelPropagation','LabelSpreading']
    data = request.get_json()['data']
    print(data)
    selected = data['anno']
    file = data['name']
    refinerID = data['refiner']

    if refinerID == 0:
        refined = EagerRefine(selected, file)
    elif refinerID == 1:
        refined = LazyRefine(selected, file)
    elif refinerID == 2:
        refined = LPARefine(selected, file)
    else:
        refined = LPARefine(selected, file, use_model=LabelSpreading)
    if len(refined) == 0:
        response = jsonify({"refined": refined, 'success': False, 'message': "Lasso is too small."})
    else:
        response = jsonify({"refined": refined, 'success': True, 'message': "Refining by {}".format(Refiners[refinerID])})
    return response


@app.route('/upload', methods=['POST'])
def upload():
    data = request.get_json()
    selected = data['selected']
    obs = data['obs']
    refined = EagerRefine(selected, obs)
    return jsonify({"refined": refined,'success': True})


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001)
