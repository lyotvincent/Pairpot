import {
  DownloadOutlined,
  EyeOutlined,
  FileTextOutlined,
  LinkOutlined,
  SearchOutlined,
} from '@ant-design/icons'
import React, { useEffect, useState, useImperativeHandle } from 'react'
import { Statistic, ConfigProvider, theme, Row, Col, Table } from 'antd'
import TextCollapse from './TextCollapse'
import './bg.scss'

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
    st: 'ST data',
    sc: 'SC data',
  },
  {
    key: '2',
    attr: 'Title',
    st: 'ST title',
    sc: 'SC title',
  },
  {
    key: '3',
    attr: 'Contributor(s)',
    st: 'ST Contributor(s)',
    sc: 'SC Contributor(s)',
  },
  {
    key: '4',
    attr: 'Summary',
    st: 'ST Summary',
    sc: 'SC Summary',
  },
  {
    key: '5',
    attr: 'Species',
    st: 'ST Species',
    sc: 'SC Species',
  },
  {
    key: '6',
    attr: 'Tissues',
    st: 'ST Tissues',
    sc: 'SC Tissues',
  },
  {
    key: '7',
    attr: 'Technologies',
    st: 'ST Technologies',
    sc: 'SC Technologies',
  },
  {
    key: '8',
    attr: 'Contacts',
    st: 'ST Contacts',
    sc: 'SC Contacts',
  },
  {
    key: '9',
    attr: 'Citation',
    st: 'ST Citation',
    sc: 'SC Citation',
  },
  {
    key: '10',
    attr: 'Accessions',
    st: 'ST Accessions',
    sc: 'SC Accessions',
  },
]

const DatasetDetails = ({onRef}) => {
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

  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": setState, // Trigger for useEffect
  }))

  var dataIndex = ["dataset_id", "title", "contributors", "summary", "species", "tissues", "technologies", "contacts", "citation", "accessions"]
  useEffect(() => {
    if (typeof state !== 'undefined') {
      let stateKeys = Object.keys(state)
      let set_sc = stateKeys.includes('sc')
      if (stateKeys.includes('st')) {
        let st_state = state.st
        let _dataSource = dataIndex.map((item, index) => {
          let _data = Data[index]
          _data.st = item === 'summary' ? <TextCollapse text={st_state[item]} /> : st_state[item]
          if (set_sc) {
            let sc_state = state.sc
            _data.sc = item === 'summary' ? <TextCollapse text={sc_state[item]} /> : sc_state[item]
          }
          return _data
        })
        setDataSource(_dataSource)
        setTitle(st_state['title'])
        // setState({
        //   Spots: st_state['spots'],
        //   Samples: st_state['n_samples'],
        //   Cells: st_state['cells'],
        //   Genes: st_state['genes'],
        // })
      }
    }
  }, [state])

  return (
    <>
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
        {/* <Row style={{ marginBottom: 20 }} justify="space-evenly">
          <Col span={4} offset={1}>
            <Statistic title="Samples" value={state.Samples} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Spots" value={state.Spots} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Genes" value={state.Genes} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Cells" value={state.Cells} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
        </Row> */}
        <Table
          columns={Columns}
          dataSource={dataSource}
          bordered
          size="small"
          pagination={false}></Table>

      </ConfigProvider>
    </>
  )
}
export default DatasetDetails