# -*- coding=utf-8 -*-
import hashlib
import re
import pandas as pd
from pandas.errors import ParserError
import logging
import traceback

from xlrd import XLRDError



def to_md5(pwd):
    encryption = hashlib.md5()
    if pwd:
        encryption.update(pwd.encode('utf-8'))
        return encryption.hexdigest()
    return False

def is_id_number(id_number):
    if len(id_number) != 18 and len(id_number) != 15:
        print('身份证号码长度错误')
        return False
    regularExpression = "(^[1-9]\\d{5}(18|19|20)\\d{2}((0[1-9])|(10|11|12))(([0-2][1-9])|10|20|30|31)\\d{3}[0-9Xx]$)|" \
                        "(^[1-9]\\d{5}\\d{2}((0[1-9])|(10|11|12))(([0-2][1-9])|10|20|30|31)\\d{3}$)"

    if re.match(regularExpression, id_number):
        if len(id_number) == 18:
            n = id_number.upper()
            # 前十七位加权因子
            var = [7, 9, 10, 5, 8, 4, 2, 1, 6, 3, 7, 9, 10, 5, 8, 4, 2]
            # 这是除以11后，可能产生的11位余数对应的验证码
            var_id = ['1', '0', 'X', '9', '8', '7', '6', '5', '4', '3', '2']

            sum = 0
            for i in range(0, 17):
                sum += int(n[i]) * var[i]
            sum %= 11
            if (var_id[sum]) != str(n[17]):
                print("身份证号规则核验失败，校验码应为", var_id[sum], "，当前校验码是：", n[17])
                return False
        return True
    else:
        return False

def read_excel(filepath, **read_param):
    """ 通过pd.read_excel 或 pd.read_csv 获取数据 """
    if filepath.endswith('csv'):
        return pd.read_csv(filepath, **read_param)
    else:
        return pd.read_excel(filepath, **read_param)
    
def upload_excel(filename, name_dict):
    """接受上传的excel并转为df
    :param str filename : 文件名
    :param dict name_dict: 对应的列名称
    :return DataFrame
    """
    try:
        df = read_excel(filename, dtype=str)
        # 去除表头空格
        col_names = df.columns.tolist()
        for index,value in enumerate(col_names):
            col_names[index]= value.strip()
        df.columns=col_names 
        if df.empty:
            raise Exception('上传文件为空，请重新上传！')
        df.rename(columns=name_dict, inplace=True)
        df = df[list(name_dict.values())].dropna(how='all', axis=0)
    except Exception as err:
        logging.error(traceback.format_exc())
        if str(err) in ["Can't find workbook in OLE2 compound document",
                        "Workbook is encrypted"]:
            raise ValueError("xls/xlsx文件已经加密,请解除加密后重新上传")
        if isinstance(err, (XLRDError, ParserError, UnicodeDecodeError, KeyError)):
            raise ValueError("导入文件格式与模板需要的文件格式不一致")
        if "Error tokenizing data" in str(err):
            raise ValueError("解析失败，可能是分隔符问题，请调整分隔符参数后重试")
        # if isinstance(err, UnicodeDecodeError):
        #     raise ValueError("解析失败，请确认该文件类型正常后,调整编码参数进行重试")
        raise err
    else:
        return df