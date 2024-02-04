# -*- coding=utf-8 -*-
"""错误类型
"""
from flask import jsonify

""" HTTP状态码配置
* _HC_SUCCESS:          访问成功
* _HC_PARAMS_ERROR:     参数异常
* _HC_STATUS_DOING:     异步执行未完成
* _HC_NODATA:           返回结果为空
* _HC_ERROR:            访问失败
* _HTTP_ERROR_API:      HTTP状态码API配置
"""
_HC_SUCCESS = 200
_HC_PARAMS_ERROR = 401
_HC_STATUS_DOING = 202
_HC_NODATA = 203
_HC_ERROR = 500
_HC_SYSTEM_ERROR = 503
_HTTP_ERROR_API = {
    _HC_PARAMS_ERROR: ('执行失败, 参数有误', None, {}),
    _HC_STATUS_DOING: ('执行中', None, {}),
    _HC_NODATA: ('执行结果为空', None, {}),
    _HC_ERROR: ('执行失败, 系统异常中断', None, {}),
    _HC_SYSTEM_ERROR: ('存在待修复问题', None, {}),
}

class BaseError(Exception):
    """自定义错误类型"""

    def __init__(self, code, description=None, model=None, **kwargs):
        self.code = code
        self.description = description if description else _HTTP_ERROR_API.get(code, '未知错误类型')
        # 返回信息
        self.result = (jsonify({'message': self.description}), self.code)


class ParamsError(BaseError):
    """ 执行失败, 参数有误 """

    def __init__(self, description=None):
        """  初始化 """
        super().__init__(_HC_PARAMS_ERROR, description)


class ExcelError(BaseError):
    """ Excel相关错误 """

    def __init__(self, description=None, **kwargs):
        """  初始化
        :param str code: 错误类型
        :param str description: 错误说明
        """
        if description == None:
            code, description = self.get_description(**kwargs)
        else:
            code = _HC_ERROR

        super().__init__(code, description)

    def get_description(self, operate='read', file_name=None, error=None, **kwargs):
        """  获取错误说明
        :param str operate: 操作方式
        :param str file_name: 文件名称
        :param str err: 错误类型
        """
        if operate == 'read':
            description = '文件' + ('[%s]' % file_name if file_name != None else '') \
                          + '读取失败'

            # 补充异常说明
            # csv: 文件无法区分各种异常类型, 报错皆为解码失败
            # xls、xlsx:
            #   1) 编码问题
            #   2) 文件加密问题
            if error != None:
                if file_name != None and file_name.endswith('csv'):
                    if "'utf-8' codec can't decode" in error:
                        description += ', 需检查: 1.编码为utf-8 2.文件未加密'
                else:
                    if "'utf-8' codec can't decode" in error:
                        description += ', 编码格式需为utf-8'
                    elif "Can't find workbook in OLE2 compound document" in error or "Workbook is encrypted" in error:
                        description += ', 无法读取加密文件, 请解除加密后重新上传'

            return _HC_PARAMS_ERROR, description
        elif operate == 'write':
            description = '文件' + ('[%s]' % file_name if file_name != None else '') \
                          + '写入失败'

            # 补充异常说明
            # 1) No such file or directory
            # 2) Permission denied
            if error != None:
                if 'No such file or directory' in error:
                    description += ', 写入目录不存在'
                elif 'Permission denied':
                    description += ', 无权限写入文件'

            return _HC_SYSTEM_ERROR, description
        else:
            return _HC_ERROR, None


class StatusError(BaseError):
    """ 执行中 """

    def __init__(self, description=None):
        """  初始化 """
        super().__init__(_HC_STATUS_DOING, description)


class NoDataError(BaseError):
    """ 执行结果为空 """

    def __init__(self, description):
        """  初始化 """
        super().__init__(_HC_NODATA, description)


class SysError(BaseError):
    """ 系统执行错误 
    例如: 无读取权限等
    """

    def __init__(self, description):
        """  初始化 """
        super().__init__(_HC_ERROR, description)
