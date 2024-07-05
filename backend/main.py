from flask import Flask, request, jsonify, send_file
import sqlite3
import time
from server.refine import *
from server.deconvolution import *
from id import builtID
app = Flask(__name__)
app.debug = True

@app.after_request
def add_cors_headers(response):
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type, x-requested-with')
    return response


@app.route('/hello', methods=['GET'])
def hello_world():
    return 'Hello, World!'


@app.route('/api/add', methods=['POST'])
def add_numbers():
    data = request.get_json()
    num1 = data['num1']
    num2 = data['num2']
    result = num1 + num2
    return jsonify({'result': result})

@app.route('/example', methods=['GET'])
def example():

    # get study id
    id = request.args.get("id")
    if id is None:
        id = 'STDS0000235'

    # connect to database
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # fetch attributes of datasets
    cursor.execute("PRAGMA table_info(datasets)")
    attributes = [row[1] for row in cursor.fetchall()]

    # get data from through a specific ID
    cursor.execute(f"SELECT * FROM datasets WHERE dataset_id = \"{id}\"")
    row1 = cursor.fetchall()

    scid = row1[0][-1]
    cursor.execute(f"SELECT * FROM datasets WHERE dataset_id = \"{scid}\"")
    row2 = cursor.fetchall()
    # return query results
    response = jsonify({'attributes': attributes, 'data': [row1[0], row2[0]]})
    
    # close connection
    cursor.close()
    conn.close()

    return response



@app.route('/datasets', methods=['GET'])
def init_dataset():

    # connect to database
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # fetch attributes of datasets
    cursor.execute("PRAGMA table_info(datasets)")
    attributes = [row[1] for row in cursor.fetchall()]

    # fetch all data
    cursor.execute("SELECT * FROM datasets")
    rows = cursor.fetchall()

    # return query results
    response = jsonify({'attributes': attributes, 'data': rows})

    # close connection
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
    selected = data['anno']
    type = data['type']
    if "id" in data.keys():
        id = data["id"]
        if id is not None:
            if id[:2] == "ST":
                id = id[-3:]
        if id not in builtID:
            id = '235'
    else:
        id='235'
    ft = "sc_sampled" if type == 'sc' else 'sp_deconv'
    file = f'./resources/{id}/{ft}.h5ad'
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
        response = jsonify({"refined": refined,
                            'success': False,
                            'message': "Lasso is too small.",
                            'endtime': int(time.time()*1000)})
    else:
        response = jsonify({"refined": refined,
                            'success': True,
                            'message': f"Refined by {Refiners[refinerID]} successfully.",
                            "endtime": int(time.time()*1000),
                            })
    return response

@app.route('/deconv', methods=['POST'])
def deconv():
    data = request.get_json()['data']
    selected = data['anno']
    if "id" in data.keys():
        id = data["id"]
        if id is not None:
            if id[:2] == "ST":
                id = id[-3:]
        if id not in builtID:
            id = '235'
    else:
        id='235'
    scfile = f'./resources/{id}/sc_sampled.h5ad'
    spfile = f'./resources/{id}/sp_deconv.h5ad'
    print(id)
    try:
        prop = NNLSDeconvolution(selected, scfile, spfile)
        response = jsonify({"props": prop,
                            'success': True,
                            'message': "Refined by PairView successfully.",
                            "endtime": int(time.time()*1000),
                            })
    except Exception as e:
        response = jsonify({
                    'success': False,
                    'message': f"Refined by PairView Failed.{e}",
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

@app.route('/query', methods=['GET'])
def query():
    id = request.args.get("id")
    if id is not None:
        if id[:2] == "ST":
            id = id[-3:]
    print(id)
    if id not in builtID:
        id = '235'
    return send_file(f'./resources/{id}/sp_meta.h5ad')

@app.route('/query/sp', methods=['GET'])
def query_sp():
    id = request.args.get("id")
    if id is not None:
        if id[:2] == "ST":
            id = id[-3:]
    print(id)
    if id not in builtID:
        id = '235'
    return send_file(f'./resources/{id}/sp_meta.h5ad')

@app.route('/query/sc', methods=['GET'])
def query_sc():
    id = request.args.get("id")
    if id is not None:
        if id[:2] == "ST":
            id = id[-3:]
    print(id)
    if id not in builtID:
        id = '235'
    return send_file(f'./resources/{id}/sc_meta.h5ad')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5522)
