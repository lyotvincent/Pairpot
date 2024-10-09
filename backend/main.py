from flask import Flask, request, jsonify, send_file
from gevent import pywsgi
import sqlite3
import time
import re
from server.refine import *
from server.deconvolution import *
from id import builtID
from sklearn.feature_extraction.text import TfidfVectorizer
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

    scid = row1[0][-2]
    cursor.execute(f"SELECT * FROM datasets WHERE dataset_id = \"{scid}\"")
    row2 = cursor.fetchall()
    # return query results
    response = jsonify({'attributes': attributes, 'data': [row1[0], row2[0]]})
    
    # close connection
    cursor.close()
    conn.close()

    return response


@app.route('/search_key', methods=['GET'])
def search_dataset():
    type = request.args.get('type')
    content = request.args.get('content')
    # print("request:",type,content)
    # connect
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # Categorical discussions
    if(type == 'all'):
        query_content = f"SELECT * FROM datasets"
        cursor.execute(query_content)
        res_list = cursor.fetchall()
        response = jsonify({'data':res_list,'type':type})
        return response

    if(type == 'key'):
        id_list = find_dataset(content)
    elif(type == 'id'):
        id_list = [content] # only one
    elif(type == 'num'):
        if(len(content)>3):
            id_list = []
        else:
            std_id = "STDS"
            scd_id = "SCDS"
            for i in range(7- len(content)):
                std_id += '0'
                scd_id += '0'
            std_id += content
            scd_id += content
            id_list = [std_id, scd_id]
    # query db
    res_list = []
    for id in id_list:
        query_content = f"SELECT * FROM datasets WHERE dataset_id = '{id}'"
        cursor.execute(query_content)
        temp = cursor.fetchall()
        if len(temp) > 0: # must exist
            res_list.append(temp[0])
    response = jsonify({'data':res_list,'type':type})
    return response


def find_dataset(key):
    key = key.replace("'", "")
    # keys = key.split(" ")
    # title > summary > overall_design > species/Tissues/Platforms/disease
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()
    query_list = ['title','summary','overall_design','species','Tissues','disease','Platforms']
    id_list = []
    for query_label in query_list:
        query_content = f"SELECT * FROM datasets WHERE {query_label} Like '%{key}%'"
        cursor.execute(query_content)
        res = cursor.fetchall()
        if res != None:
            for item in res:
                id_list.append(item[1]) # record the id
    # There is no result 
    if len(id_list) == 0:
        return [] # TODO:Fuzzy queries
    # count the frequency
    count_dict = {}
    for item in id_list:
        if item in count_dict:
            count_dict[item] += 1
        else:
            count_dict[item] = 1
    # find the max value
    max_count = max(count_dict.values())
    # find items with the max value
    res_list = []
    unique_list = list(set(id_list))
    for item in unique_list:
        if count_dict[item] == max_count:
            res_list.append(item)
    return res_list  


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

@app.route('/global_words', methods=['GET'])
def get_global_wordcloud():
    # connect to database
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # query content
    query_pair = {
        'datasets': ["title", "species",	"tissues", "organ_parts", "cell_types", "summary", "overall_design"],
        'pair_qualities': ["quality", "description"],
        'sc_datasets': ["title", "species",	"tissues", "organ_parts", "cell_types", "summary", "overall_design"],
        'species_table': ["species_name", "simple_name"],
        'technology_table': ["point_address", "spot_configuration", "strategy"],
        'tissue_table': ["tissue_name"]
    }

    res = []
    for table_name, table_attr in query_pair.items():
        for i in range(len(table_attr)):
            temp_query = f"SELECT {table_attr[i]} FROM {table_name}"
            cursor.execute(temp_query)
            res.append(cursor.fetchall())

    text = []
    for i in range(len(res)):
        for j in range(len(res[i])):
            item_str = str(res[i][j][0])
            item_str = re.sub(r'[^\w\s]', ' ', item_str)
            words = item_str.split(" ")
            text.append(words)
    
    corpus = [' '.join(sentence) for sentence in text]
    # print(corpus)

    custom_stop_words = [
        'none','we','were',
        'nan','homo','sapiens','mus','musculus',
        'and', 'if', 'or', 'the', 'to', 'in', 'of', 'a', 'an', 'is', 'was', 'for', 'on', 'with', 'at', 'by', 'from', 'as', 'that', 'which', 'it', 'are', 'has', 'have', 'had', 'be', 'will', 'would', 'can', 'could', 'may', 'might', 'should'
    ]

    # init TfidfVectorizer
    vectorizer = TfidfVectorizer(stop_words=custom_stop_words)
    # TF-IDF Matrix
    X = vectorizer.fit_transform(corpus)
    # TF-IDF weghts
    feature_names = vectorizer.get_feature_names_out()
    dense = X.todense()
    denselist = dense.tolist()
    df = pd.DataFrame(denselist, columns=feature_names)
    word_tfidf = df.sum(axis=0)
    #  to dict
    word_dict = word_tfidf.to_dict()

    return word_dict


@app.route('/filted_words', methods=['POST'])
def get_filted_wordcloud():
    response = request.get_json()['params']
    print(response)
    # save : title, summary, overall_design, species,organ_parts, cell_types
    # res = [[item[2], item[3], item[4], item[5], item[6],item[12], item[16], item[17]] for item in res]
     # changed param
    label = None if "label" not in response.keys() else response["label"]
    item =  None if "item" not in response.keys() else response["item"]
    src = [] if "src" not in response.keys() else response["src"]
    print(label, item, src)

    full_content = query_content = "SELECT title, summary, overall_design, species,organ_parts, cell_types FROM datasets"
    
    # step1: connect
    conn = sqlite3.connect('resources/STpair.db')
    cursor = conn.cursor()

    # step2: query
    if label == "all":  # using label & item
        query_content = full_content
        cursor.execute(query_content)
    elif label is not None and item is not None:
        query_content = f"SELECT title, summary, overall_design, species,organ_parts, cell_types FROM datasets where {label} like '%{item}%'"
        cursor.execute(query_content)
    elif len(src) > 0:  # using src 
        cursor.execute(f"SELECT COUNT(*) FROM datasets")  # fetch the data length
        row_count = cursor.fetchone()[0]
        if len(src) >= row_count:  # return all if src array is too large
            query_content = full_content
            cursor.execute(query_content)
        else:
            query_content = f"SELECT title, summary, overall_design, species,organ_parts, cell_types FROM datasets where dataset_id IN ({','.join(['?' for _ in src])})"
            cursor.execute(query_content, src)
    else:  # return message error
        cursor.close()
        conn.close()
        return jsonify({"error": "Missing required labels or filter words"}), 400
    res = cursor.fetchall()

    # step3: calculate TF-IDF    
    text = []
    for i in range(len(res)):
        item_str = str(res[i])
        item_str = re.sub(r'[^\w\s]', ' ', item_str)
        words = item_str.split(" ")
        text.append(words)
        # print(i)
    # print(text)
    corpus = [' '.join(sentence) for sentence in text]
    # print(corpus)
    custom_stop_words = [
        'none','we','were',
        'nan','homo','sapiens','mus','musculus',
        'and', 'if', 'or', 'the', 'to', 'in', 'of', 'a', 'an', 'is', 'was', 'for', 'on', 'with', 'at', 'by', 'from', 'as', 'that', 'which', 'it', 'are', 'has', 'have', 'had', 'be', 'will', 'would', 'can', 'could', 'may', 'might', 'should'
    ]
    # init TfidfVectorizer
    vectorizer = TfidfVectorizer(stop_words=custom_stop_words)
    # TF-IDF Matrix
    X = vectorizer.fit_transform(corpus)
    # TF-IDF weghts
    feature_names = vectorizer.get_feature_names_out()
    dense = X.todense()
    denselist = dense.tolist()
    df = pd.DataFrame(denselist, columns=feature_names)
    word_tfidf = df.sum(axis=0)
    #  to dict
    word_dict = word_tfidf.to_dict()
    # print(word_dict)

    return jsonify({"data": word_dict})


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

@app.route('/query/spgz', methods=['GET'])
def query_spgz():
    id = request.args.get("id")
    if id is not None:
        if id[:2] == "ST":
            id = id[-3:]
    print(id)
    if id not in builtID:
        id = '235'
    return send_file(f'./resources/{id}/sp_meta.h5ad.zip')

@app.route('/query/scgz', methods=['GET'])
def query_scgz():
    id = request.args.get("id")
    if id is not None:
        if id[:2] == "ST":
            id = id[-3:]
    print(id)
    if id not in builtID:
        id = '235'
    return send_file(f'./resources/{id}/sc_meta.h5ad.zip')

if __name__ == '__main__':
    server = pywsgi.WSGIServer(('127.0.0.1', 5522), app)
    server.serve_forever()
    # get_global_wordcloud()
