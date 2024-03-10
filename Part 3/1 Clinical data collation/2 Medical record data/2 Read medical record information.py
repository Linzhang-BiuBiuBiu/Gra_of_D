# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 10:12:05 2021

@author: zl
"""

import os
from pathlib import Path
import pandas as pd
import re
import docx2txt

struct = {"姓名": "value","入院日期": "value","性别": "value","出院日期": "value","年龄": "value","住院天数": "value","入院情况": "value",
    "查体": "value",
    "mmHg": "value",
    "BP": "value",
    "NP": "value",
    "EF%": "value",
    "入院诊断": "value",
    "诊疗经过": "value",
    "心脏超声":"value",
    "心脏彩超":"value",
    "BNP": "value",
    "EF": "value",
    "冠状动脉造影":"value",
    "出院诊断": "value",
    "出院情况":"value",
    "出院医嘱": "value",
}


# 构建函数
def read_key_value(obj, buff):
    key_list = []
    for k in obj.keys():
        key_list.append(k)
    key_list.append("$")
    for n in range(len(key_list) - 1):
        key = key_list[n]
        # key_next = key_list[n+1]
        key_with_space = "\s*".join(list(key))
        # 判断key是否存在，若不存在，value直接设置为空值
        key_current = re.findall(key_with_space, buff)
        if len(key_current) == 0:
            obj[key] = ""
            continue
        # 若key存在，寻找下一个key作为终止符
        for next_n in range(n + 1, len(key_list)):
            next_key = key_list[next_n]
            next_key_with_space = "\s*".join(list(next_key))
            next_key_current = re.findall(next_key_with_space, buff)
            # 若下一个值不是终止符，也没有存在病例中的话，寻找下一个key
            if len(next_key_current) == 0 and next_key != "$":
                continue
            # 找到key与next_key，把正则表达式的值给value
            pat = re.compile(key_with_space + "[:|：]{0,1}(.*?)" + next_key_with_space, re.S)
            info_list = re.findall(pat, buff)
            content = "".join(info_list)
            # 判断key是否包含层级结构，如果包含，或的内容再次使用key进行拆解
            if isinstance(obj[key], dict):
                read_key_value(obj[key], content)
            else:
                obj[key] = content
            break
    return None


# 数据转换成DataFame
cols = ["文件名","入院日期","性别","出院日期","年龄","住院天数","入院情况","查体","mmHg","BP","NP","EF%","入院诊断","诊疗经过","BNP","EF","心脏超声","心脏彩超","冠状动脉造影","出院诊断","出院情况","出院医嘱"]

df = pd.DataFrame(columns=cols)
file_path = Path(r"D:\onedrive\桌面\Clinical_data_summary\入院+出院+手术-v2\出院")

# 用rglob获取docx文件
for file in file_path.rglob('*.docx'):
    # 读取每个病历
    content = docx2txt.process(file)
    print(file)
    # 调用read_key_value函数，依次读取key值和value值
    read_key_value(struct, content)
    record = {}
    for c in cols:
        if c == '文件名':
            record[c] = file.name
        else:
            record[c] = struct[c]
    # print(record)
    df = df.append(record, ignore_index=True)



df.to_excel(r"D:\onedrive\桌面\Clinical_data_summary\入院+出院+手术-v2\住院.xlsx")
