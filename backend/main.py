from flask import Flask, request, jsonify
import sqlite3
import time
from server.refine import *
app = Flask(__name__)
app.debug = True

@app.after_request
def add_cors_headers(response):
    response.headers.add('Access-Control-Allow-Origin', 'http://localhost:6634')
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
    starttime = data['starttime']

    if refinerID == 0:
        refined = EagerRefine(selected, file)
    elif refinerID == 1:
        refined = LazyRefine(selected, file)
    elif refinerID == 2:
        refined = LPARefine(selected, file)
    else:
        refined = LPARefine(selected, file, use_model=LabelSpreading)
    if len(refined) == 0:
        response = jsonify({"refined": refined,
                            'success': False,
                            'message': "Lasso is too small.",
                            'endtime': int(time.time()*1000)})
    else:
        response = jsonify({"refined": refined,
                            'success': True,
                            'message': "Refined by {} successfully.".format(Refiners[refinerID]),
                            "endtime": int(time.time()*1000),
                            })
    return response


@app.route('/upload', methods=['POST'])
def upload():
    data = request.get_json()
    selected = data['selected']
    obs = data['obs']
    refined = EagerRefine(selected, obs)
    return jsonify({"refined": refined,'success': True})

@app.route('/submit', methods=['POST'])
def submit():
    info = request.get_json()
    data = info['data']
    pair_state = data.pop('pair-state', None)
    print(pair_state)
    table_name = 'datasets'

    # open the database
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    if pair_state == 'From original dataset':
        try:
            # update data
            update_query = f"UPDATE {table_name} SET has_paired = ? WHERE dataset_id = ?"
            cursor.execute(update_query, (data['has_paired'], data['has_paired']))
            conn.commit()
            message = 'Submit success.'
            state = 'success'
        except sqlite3.Error as e:
            message = print(f"Error occurs: {e}.")
            print(message)
            state = 'error'
    elif pair_state == 'From existing dataset':
        state = 'success'
        message = 'Submit success.'
    else:

        # find current id
        query = f"SELECT MAX(id) FROM {table_name}"
        cursor.execute(query)
        curr_id = cursor.fetchone()[0]
        data['id'] = curr_id + 1
        data['dataset_id'] = f"STDSA000{curr_id + 1}"

        # split species, tissues, and technologies
        for key in ['species', 'tissues', 'technologies']:
            data[key] = ';'.join(data[key])

        try:
            # insert data
            insert_query = f"INSERT INTO {table_name} ({','.join([col for col in data.keys()])}) VALUES ({','.join(['?' for _ in data.keys()])})"
            cursor.execute(insert_query, list(data.values()))
            conn.commit()
            message = 'Submit success.'
            state = 'success'
        except sqlite3.Error as e:
            message = print(f"Error occurs: {e}.")
            print(message)
            state = 'error'
    cursor.close()
    conn.close()

    return jsonify({'state': state, 'message': message})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5522)
