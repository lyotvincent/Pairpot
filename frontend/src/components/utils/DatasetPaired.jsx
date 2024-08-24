import {
  DownloadOutlined,
  EyeOutlined,
  FileTextOutlined,
  LinkOutlined,
  SearchOutlined,
} from '@ant-design/icons'
import React, { useEffect, useState, useImperativeHandle } from 'react'
import { Statistic, ConfigProvider, theme, Row, Col, Table } from 'antd'
import { Descriptions, Modal, Space, Tooltip } from 'antd'

import TextCollapse from './TextCollapse'
import './bg.scss'


const DatasetPaired = ({ descCol, descInfo, text, placement, scInfo }) => {
  const [modal1Open, setModal1Open] = useState(false)
  // console.log(scInfo)

  const Columns = [
    {
      title: 'Attributes',
      dataIndex: 'attr',
      key: 'attr',
      width: '10%',
      align: 'right',
    },
    {
      title: 'Spatial Transcriptomics data',
      dataIndex: 'st',
      key: 'st',
      width: '45%',
      align: 'center',
    },
    {
      title: 'Single-Cell data',
      dataIndex: 'sc',
      key: 'sc',
      align: 'center',
    },
  ]
  const Data = [
    {
      key: '1',
      attr: 'Dataset ID',
      st: descInfo[1],
      // sc: scInfo[1],
      sc: scInfo == null? 'SC data': scInfo[1]
    },
    {
      key: '2',
      attr: 'Title',
      st: descInfo[2],
      sc: scInfo == null? 'SC title': scInfo[2],
    },
    {
      key: '3',
      attr: 'Contributor(s)',
      st: descInfo[20],
      sc: scInfo == null? 'SC Contributor(s)': scInfo[20],
    },
    {
      key: '4',
      attr: 'Summary',
      // st: descInfo[16],
      // sc: scInfo[16],
      st: <TextCollapse text={descInfo[16]} />,
      // sc: <TextCollapse text={scInfo == null? 'SC Summary': scInfo[16]} />,
      sc: scInfo == null?  'SC Summary' : <TextCollapse text={scInfo[16]}/>

    },
    {
      key: '5',
      attr: 'Species',
      st: descInfo[3],
      sc: scInfo == null? 'SC Species': scInfo[3],
    },
    {
      key: '6',
      attr: 'Tissues',
      st: descInfo[4],
      sc: scInfo == null? 'SC Tissues': scInfo[4],
    },
    {
      key: '7',
      attr: 'Technologies',
      st: descInfo[12],
      sc: scInfo == null? 'SC Technologies': scInfo[12],
    },
    {
      key: '8',
      attr: 'Contacts',
      st: descInfo[21],
      sc: scInfo == null? 'SC Contacts': scInfo[21],
    },
    {
      key: '9',
      attr: 'Citation',
      st: descInfo[22],
      sc: scInfo == null? 'SC Citation': scInfo[22],
    },
    {
      key: '10',
      attr: 'Accessions',
      st: descInfo[23],
      sc: scInfo == null? 'SC Accessions': scInfo[23],
    },
  ]

  const [title, setTitle] = useState()
  const [state, setState] = useState()
  const [dataSource, setDataSource] = useState(Data)
  const { token } = theme.useToken()
  const formItemLayout = {
    labelCol: {
      span: 3,
    },
    wrapperCol: {
      span: 9,
    },
  }

  // useImperativeHandle(onRef, () => ({  // explode trigger for parent components
  //   "Trigger": setState, // Trigger for useEffect
  // }))

  // var dataIndex = ["dataset_id", "title", "contributors", "summary", "species", "tissues", "technologies", "contacts", "citation", "accessions"]
  // useEffect(() => {
  //   if (typeof state !== 'undefined') {
  //     let stateKeys = Object.keys(state)
  //     let set_sc = stateKeys.includes('sc')
  //     if (stateKeys.includes('st')) {
  //       let st_state = state.st
  //       let _dataSource = dataIndex.map((item, index) => {
  //         let _data = Data[index]
  //         _data.st = item === 'summary' ? <TextCollapse text={st_state[item]} /> : st_state[item]
  //         if (set_sc) {
  //           let sc_state = state.sc
  //           _data.sc = item === 'summary' ? <TextCollapse text={sc_state[item]} /> : sc_state[item]
  //         }
  //         return _data
  //       })
  //       setDataSource(_dataSource)
  //       setTitle(st_state['title'])
  //     }
  //   }
  // }, [state])

  return (
    <div>
      <Tooltip
        placement={placement}
        title="Show Paired Datasets"
        arrow={false}
        align={'center'}>
        <Space
          onClick={() => setModal1Open(true)}
          style={{ cursor: 'pointer' }}>
          <LinkOutlined />
          {text}
        </Space>
      </Tooltip>
      <Modal
        title="Paired Datasets"
        open={modal1Open}
        onCancel={() => setModal1Open(false)}
        footer={null}
        width={1200}
      >
        <ConfigProvider
          theme={{
            components: {
              Form: {
                itemMarginBottom: '0px',
              },
              Button: {
                borderColorDisabled: 'gray',
                colorTextDisabled: 'gray',
                colorBgContainerDisabled: 'rgba(0.5, 0.5, 0.5, 0.1)',
              },
              Upload: {
                padding: 0,
              },
              Statistic: {
                titleFontSize: 18,
                contentFontSize: 36,

              }
            },
          }}>
          <h3 style={{ color: 'black', marginBottom: 20 }}>{title}</h3>
          <Table
            columns={Columns}
            dataSource={dataSource}
            bordered
            size="small"
            pagination={false}></Table>

        </ConfigProvider>
      </Modal>
    </div>
  )
}
export default DatasetPaired