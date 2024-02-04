# -*- coding:utf-8 -*-
"""系统运行数据库连接工具
使用sqlalchemy的数据库相关的代码，仅包括数据库连接管理, ORM模型定义在model目录.

当前支持Oracle与Mysql两种方式

"""
import datetime
import os
from string import Template

# sqlite不支持big int autoincrement
import sqlalchemy
from sqlalchemy import (Boolean, Column, Date, DateTime, Float, ForeignKey, ForeignKeyConstraint, Index, Integer,
                        Integer as BigInteger, PrimaryKeyConstraint, Sequence, String, Table, Text,
                        UnicodeText, UniqueConstraint, create_engine, func)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.inspection import inspect
from sqlalchemy.orm import backref, mapper, relationship, scoped_session, sessionmaker
from sqlalchemy.pool import NullPool
from sqlalchemy import Text as CLOB


os.environ['NLS_LANG'] = 'SIMPLIFIED CHINESE_CHINA.UTF8'

func_to_char = func.date_format
func_to_date = func.str_to_date

#####################
# ORM 对象基本申明类#
#####################
# 函数更新了 做如下修改
# Base = declarative_base()
Base = sqlalchemy.orm.declarative_base()
# 配置项的值是一个 URL 的形式，‘mysql://user:password@ip:port/database’ ，分别是使用的数据库，登录用户，密码，ip地址，端口，数据库名。
# 例如 app.config['SQLALCHEMY_DATABASE_URI'] = 'mysql://admin:Mysql!123@127.0.0.1:3306/MyDB_one'
#SQLALCHEMY_DATABASE_URL = u"mysql+pymysql://root:123456@localhost:3306/test"
SQLALCHEMY_DATABASE_URL = u"mysql+pymysql://root:926205605518@localhost:3306/hospital_update"
# HOSTNAME = '127.0.0.1'
# PORT = '3306'
# DATABASE = 'hospital'
# USERNAME = 'root'
# PASSWORD = 'NANkai@101'
# DB_URL = 'mysql+pycharm://{}:{}@{}:{}/{}?charset=utf8'.format(USERNAME,PASSWORD,HOSTNAME,PORT,DATABASE)
# SQLALCHEMY_DATABASE_URI = DB_URL


def to_dict(self, convert=True):
    """Extend Ablity To Dict Orm Object Data"""
    attrlist = [a for a in list(self.__dict__.keys()) if not a.startswith('_')]
    data = {}
    for name in attrlist:
        attr_data = getattr(self, name, None)
        if not convert:
            data[name] = attr_data
            continue

        if isinstance(attr_data, datetime.datetime):
            attr_data = attr_data.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(attr_data, datetime.date):
            attr_data = attr_data.strftime('%Y-%m-%d')
        elif isinstance(attr_data, datetime.time):
            attr_data = attr_data.strftime('%H:%M:%S')

        data[name] = attr_data

    return data


setattr(Base, 'to_dict', to_dict)


class Database(object):
    """ Database Manager Object"""

    def __init__(self, configure, name, echo=False, pool_size=10, pool_recycle=1800,
                 poolclass=None, thread=True):
        """DB Connection Init

        :param dict configure: dblink config dict, exp:{DBNAME:DBURL}
        :param str name:dbname in configure need be connect
        :param bool echo: echo sql detail
        :param int pool_size: connection pool size
        :param int pool_recycle: recycle connection time
        :param poolclass: SQLAlchemy Pool Class, Default None
        :param bool thread: if Shared Connection between theads.Using sqlalchemy scoped_session
        """
        self.configure = configure
        extend_args = {'pool_size': pool_size}
        if poolclass == NullPool:
            extend_args = {"poolclass": NullPool}
        self.engine = create_engine(self.get_url(name), echo=echo, pool_recycle=pool_recycle, pool_pre_ping=True,
                                    encoding="utf-8", convert_unicode=False, **extend_args)
        if thread:
            self.Session = scoped_session(sessionmaker(bind=self.engine, autocommit=False,
                                                       autoflush=False))
        else:
            self.Session = sessionmaker(bind=self.engine, autocommit=False, autoflush=False)

    def __getattr__(self, name):
        if name == "session":
            return self.Session()

    def get_url(self, config_name):
        """Get DB Url From Config Dict By Name"""
        self.url = self.configure.get(config_name)
        return self.url


_SYSTEM_DB = Database({'sysdb_url': SQLALCHEMY_DATABASE_URL}, 'sysdb_url')


def simple_session(url=SQLALCHEMY_DATABASE_URL, process=False):
    if process:
        return Database({'url': url}, 'url', pool_size=1, thread=False).session
    else:
        return _SYSTEM_DB.session


def scope_session(*args, **kwargs):
    return _SYSTEM_DB.Session


def get_session(subject):
    """Get Session by inspect In Connection"""
    stat = inspect(subject)
    return stat.session


if __name__ == '__main__':
    session = simple_session()
    print(dir(session.bind))
    print(session.bind.url)
