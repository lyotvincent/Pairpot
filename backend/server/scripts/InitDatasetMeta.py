# -*- coding:utf-8 -*-
import logging
from faker import Faker
import random


fake = Faker(locale=['zh_CN'])


def init_patient(session):
    """ 初始化病人信息
    :param session:
    """

    patientList = []
    for i in range(100):
        patient = Patient(**{
            "id": random.randint(1000, 3000),
            "name": fake.name(),
            # "phone" : int(fake.phone_number()),
            "age": random.randint(8, 80),
            "fbDate": fake.date(pattern="%Y-%m-%d", end_datetime=None),
            "ryDate": fake.date(pattern="%Y-%m-%d", end_datetime=None),
            "site": random.choice(["右乳", "左乳", "双乳"])
        })

        patientList.append(patient)

    session.add_all(patientList)
    session.commit()