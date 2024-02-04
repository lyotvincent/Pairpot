# -*- coding=utf-8 -*-
import os
import logging
import traceback
import collections
import datetime
import json
import pkgutil
import importlib
from decimal import Decimal

from flask import Flask, request, g, jsonify, Blueprint
from json import JSONEncoder as BaseJSONEncoder
from flask_cors import CORS
import flask_restplus
from sqlalchemy.ext.declarative import DeclarativeMeta
from flask_jwt_extended import JWTManager, verify_jwt_in_request, get_jwt_identity

from server.database.sqlal import simple_session
from server.tools.error import BaseError


def convert(data):
    if isinstance(data, str):
        return str(data)
    elif isinstance(data, collections.Mapping):
        return dict(list(map(convert, iter(data.items()))))
    elif isinstance(data, collections.Iterable):
        return type(data)(list(map(convert, data)))
    else:
        return data


class JSONEncoder(BaseJSONEncoder):
    """Custom :class:`JSONEncoder` which respects objects that include the
    :class:`JsonSerializer` mixin.
    """

    def default(self, obj):
        if isinstance(obj, datetime.datetime):
            return obj.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(obj, datetime.date):
            return obj.strftime('%Y-%m-%d')
        elif isinstance(obj, datetime.time):
            return obj.strftime('%H:%M:%S')
        elif hasattr(obj, 'asdict') and isinstance(getattr(obj, 'asdict'), collections.Callable):
            return obj.asdict()
        # elif isinstance(obj, collections.Mapping):
        #    return dict(map(convert, obj.iteritems()))
        elif isinstance(obj, Decimal):
            return str(obj)
        elif isinstance(obj, collections.Iterable):
            return type(obj)(list(map(convert, obj)))
        elif isinstance(obj.__class__, DeclarativeMeta):
            fields = {}
            for field in [x for x in dir(obj) if
                          not x.startswith('_') and x not in ['metadata', 'query', 'query_class']]:
                data = obj.__getattribute__(field)
                try:
                    # this will fail on non-encodable values, like other classes
                    json.dumps(data)
                    fields[field] = data
                except TypeError:
                    if isinstance(data, datetime.datetime):
                        fields[field] = data.strftime('%Y-%m-%dT%H:%M:%S')
                    elif isinstance(data, datetime.date):
                        try:
                            fields[field] = data.strftime('%Y-%m-%d')
                        except ValueError:
                            fields[field] = ''

                    elif isinstance(data, datetime.time):
                        fields[field] = data.strftime('%H:%M:%S')
                    else:
                        fields[field] = None
            # a json-encodable dict
            return fields
        else:
            return super(JSONEncoder, self).default(obj)


def configure_namespace(api, package_name, package_path):
    """Register all flask_restplus Namespace instance on the specified Flask application found
    in all modules for the specified package.

    :param api: the Flask_RestPlus API
    :param package_name: the package name
    :param package_path: the package path
    """
    if not isinstance(api, flask_restplus.Api):
        raise ValueError("Must Be Flask_RESTPlus Api Object")

    all_namespace = []
    for _, name, _ in pkgutil.iter_modules(package_path):
        module = importlib.import_module('%s.%s' % (package_name, name))
        for item in dir(module):
            item = getattr(module, item)
            if isinstance(item, flask_restplus.Namespace):
                api.add_namespace(item)
            all_namespace.append(item)
    return all_namespace


def configure_blueprints(app, package_name, package_path):
    """Register all Blueprint instances on the specified Flask application found
    in all modules for the specified package.
    :param app: the Flask application
    :param package_name: the package name
    :param package_path: the package path
    """
    all_blueprints = []
    for _, name, _ in pkgutil.iter_modules(package_path):
        module = importlib.import_module('%s.%s' % (package_name, name))
        for item in dir(module):
            item = getattr(module, item)
            if isinstance(item, Blueprint):
                app.register_blueprint(item)
            all_blueprints.append(item)
    return all_blueprints


def configure_before_handlers(app):
    """Flask 服务请求前后处理"""

    @app.before_request
    def before_request():
        g.db_session = simple_session()
        # g.client_ip = request.remote_addr
        # g.message_param = {}

        if request.method == "OPTIONS":
            return

        if request.path == "/api/login":
            g.message = '登录@尝试登录'
        elif request.path == '/api/info/add':
            g.massage = '注册'
        else:
            try:
                # if request.path == "/User/RefreshToken":
                #
                # else:
                verify_jwt_in_request()
                g.user = get_jwt_identity()
                if not bool(g.user):
                    raise Exception('认证已失效!')
            except Exception as e:
                logging.error(traceback.format_exc())
                return jsonify({"error": "用户鉴权失败！"}), 401

    @app.after_request
    def after_request(response):
        """请求后执行操作审计与数据持久化 """
        try:
            g.db_session.commit()
        except Exception as e:
            logging.error(traceback.format_exc())
            g.db_session.rollback()
        finally:
            # 1. 判断是否为get_parser.parse_args未找到较好的判断条件
            # 2. 未能在errorhandler中捕获到400异常
            if response.status_code == 400:
                info = response.get_json()
                if 'message' in info and 'errors' in info:
                    param_names = []
                    for error in info['errors'].values():
                        param_names.append(error.split(' ')[0])
                    response = jsonify({'message': '参数校验失败, 请检查参数[%s]' % (', '.join(param_names))})
                    response.status_code = 400

        response.headers.add('Access-Control-Allow-Credentials', 'true')
        response.headers.add('Access-Control-Allow-Origin', request.headers.get("Origin"))
        response.headers.add('Access-Control-Allow-Headers',
                             'Content-Type,Authorization,x-session-token,X-Requested-With')
        response.headers.add('Access-Control-Allow-Methods', 'HEAD, OPTIONS, GET,PUT,POST,DELETE')

        return response

    @app.teardown_request
    def teardown_request(exception):
        g.db_session = simple_session()
        if g.db_session:
            g.db_session.close()

    @app.errorhandler(Exception)
    def APIException(error):
        if isinstance(error, BaseError):
            # 当前不同类型无特殊处理
            return error.result
        elif isinstance(error, (ValueError, RuntimeError)):
            return jsonify({'message': str(error)}), 500
        else:
            logging.error(traceback.format_exc())
            return jsonify({'message': '请联系系统管理员'}), 500


def configure_jwtmanager(app):
    # token生效时间(单位:min)
    _TOKEN_ACCESS_TIME = 700
    # # token刷新时限
    # _TOKEN_REFRESH_TIME = 1440
    # token签名的加密方式
    _ALGORITHM = os.environ.get('_ALGORITHM', 'HS256')

    # TODO 上线修改
    # 对称签名算法的密钥
    _TOKEN_SECRET = 'jjdw'
    # 非对称算法的公钥
    _TOKEN_PUB_SECRET = ''
    # 非对称算法的私钥
    _TOKEN_PER_SECRET = ''

    # 用户声明是否应包含在刷新令牌中
    _CLAIMS_IN_REFRESH = False
    # token在headers中的key
    _TOKEN_HEADER_NAME = 'x-session-token'
    # token在headers中的value前缀
    _JWT_HEADER_TYPE = ''
    # token版本号
    _TOKEN_VERSION = 1.1

    app.config['JWT_ACCESS_TOKEN_EXPIRES'] = datetime.timedelta(minutes=_TOKEN_ACCESS_TIME)
    # app.config['JWT_REFRESH_TOKEN_EXPIRES'] = datetime.timedelta(minutes=_TOKEN_REFRESH_TIME)
    app.config['JWT_ALGORITHM'] = _ALGORITHM
    if app.config['JWT_ALGORITHM'] == 'HS256':
        app.config['JWT_SECRET_KEY'] = _TOKEN_SECRET
    else:
        app.config['JWT_PUBLIC_KEY'] = _TOKEN_PUB_SECRET
    app.config['JWT_PRIVATE_KEY'] = _TOKEN_PER_SECRET
    app.config['JWT_CLAIMS_IN_REFRESH_TOKEN'] = _CLAIMS_IN_REFRESH
    app.config['JWT_HEADER_NAME'] = _TOKEN_HEADER_NAME
    app.config['JWT_HEADER_TYPE'] = _JWT_HEADER_TYPE

    jwt = JWTManager(app)


def create_app(package_name, package_path):
    """Create Flask Application

    :param str package_name:
    :param str package_path:
    """
    app = Flask(__name__)

    # from .settings import Config
    # app.config.from_object(config)
    #
    # configure_blueprints(app, package_name, package_path)
    #
    # app.wsgi_app = HTTPMethodOverrideMiddleware(app.wsgi_app)

    CORS(app)
    app.json_encoder = JSONEncoder

    # before_request、after_request..
    configure_before_handlers(app)
    # setting jwt
    configure_jwtmanager(app)

    return app


app = create_app(__name__, __path__)
environ = {'wsgi.version': (1, 0), 'wsgi.input': '', 'REQUEST_METHOD': 'GET', 'PATH_INFO': '/',
           'SERVER_NAME': 'itcast server', 'wsgi.url_scheme': 'http', 'SERVER_PORT': '80'}  # 模拟解析客户端请求之后的wsgi字典数据

with app.request_context(environ):
    api = flask_restplus.Api(app, "YINSHO", title="XX系统", description="系统后端")
    configure_namespace(api, "rptserver.views", [os.path.dirname(__file__)])
