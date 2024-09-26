import {
  DownloadOutlined,
  EyeOutlined,
  FileTextOutlined,
  LinkOutlined,
  SearchOutlined,
} from '@ant-design/icons'
import React, { useEffect, useState, useImperativeHandle } from 'react'
import { Statistic, ConfigProvider, theme, Row, Col, Table, Space, Button, Input, Select, Radio } from 'antd'
import TextCollapse from './TextCollapse'
import './bg.scss'
const { TextArea } = Input

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
  {
    key: '11',
    attr: 'Processed files',
    st: 'ST processed files',
    sc: 'SC processed files',
  },
]

const DatasetDetails = ({ onRef }) => {
  const [title, setTitle] = useState()
  const [state, setState] = useState()
  const [stid, setStid] = useState()
  const [bulkText, setBulkText] = useState()
  const [bulkTool, setBulkTool] = useState("wget")
  const [statistic, setStatistic] = useState({
    Spots: 0,
    Slices: 0,
    Cells: 0,
    Genes: 0,
  })
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

  const BulkDownload = (tool, stid) => {
    if (tool === 'curl') {
      return (
        `# Processed Files\ncurl -o http://src.bioxai.cn/${stid}/sp_deconv.h5ad\ncurl -o http://src.bioxai.cn/${stid}/sp_meta.h5ad.zip\ncurl -o http://src.bioxai.cn/${stid}/sc_sampled.h5ad\ncurl -o http://src.bioxai.cn/${stid}/sc_meta.h5ad.zip`
      )
    } else {
      return (
        `# Processed Files\nwget http://src.bioxai.cn/${stid}/sp_deconv.h5ad\nwget http://src.bioxai.cn/${stid}/sp_meta.h5ad.zip\nwget http://src.bioxai.cn/${stid}/sc_sampled.h5ad\nwget http://src.bioxai.cn/${stid}/sc_meta.h5ad.zip`
      )
    }
  }

  var dataIndex = ["dataset_id", "title", "contributors", "summary", "species", "tissues", "technologies", "contacts", "citation", "accessions"]
  useEffect(() => {
    console.log(state)
    if (typeof state !== 'undefined') {

      let stateKeys = Object.keys(state)
      let set_sc = stateKeys.includes('sc')
      if (stateKeys.includes('st')) {
        let st_state = state.st
        let _stat = {
          Spots: st_state['spots'],
          Slices: st_state['n_samples'],
          Genes: st_state['genes'],
        }
        // set processed files
        let _procFiles = Data[Data.length - 1]
        _procFiles.st = <Space>
          <Button href={`http://src.bioxai.cn/${st_state['dataset_id'].slice(-3)}/sp_deconv.h5ad`} style={{ textDecorationLine: 'none' }}>sp_deconv.h5ad</Button>
          <Button href={`http://src.bioxai.cn/${st_state['dataset_id'].slice(-3)}/sp_meta.h5ad.zip`} style={{ textDecorationLine: 'none' }}>sp_meta.h5ad.zip</Button>
        </Space>
        if (set_sc) {
          _procFiles.sc = <Space>
            <Button href={`http://src.bioxai.cn/${st_state['dataset_id'].slice(-3)}/sc_sampled.h5ad`} style={{ textDecorationLine: 'none' }}>sc_sampled.h5ad</Button>
            <Button href={`http://src.bioxai.cn/${st_state['dataset_id'].slice(-3)}/sc_meta.h5ad.zip`} style={{ textDecorationLine: 'none' }}>sc_meta.h5ad.zip</Button >
          </Space>
        }
        // set bulk download

        let _dataSource = dataIndex.map((item, index) => {
          let _data = Data[index]
          _data.st = item === 'summary' ? <TextCollapse text={st_state[item]} /> : st_state[item]
          if (set_sc) {
            let sc_state = state.sc
            _stat = {
              ..._stat,
              Cells: sc_state['cells'],
            }
            _data.sc = item === 'summary' ? <TextCollapse text={sc_state[item]} /> : sc_state[item]
          }
          return _data
        })
        _dataSource.push(_procFiles)
        setStatistic(_stat)
        setDataSource(_dataSource)
        setBulkText(BulkDownload(bulkTool, st_state['dataset_id'].slice(-3)))
        setTitle(st_state['title'])
        setStid(st_state['dataset_id'].slice(-3))
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
        <Row style={{ marginBottom: 20 }} justify="space-evenly">
          <Col span={4} offset={1}>
            <Statistic title="Slices" value={statistic.Slices} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Spots" value={statistic.Spots} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Genes" value={statistic.Genes} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
          <Col span={4} offset={1}>
            <Statistic title="Cells" value={statistic.Cells} valueStyle={{ color: token.colorPrimaryActive }} />
          </Col>
        </Row>
        <Table
          columns={Columns}
          dataSource={dataSource}
          bordered
          size="small"
          pagination={false}>
        </Table>
        <Space style={{ marginTop: 10 }}>
          <h5 style={{ color: 'black', marginTop: 10 }}>Bulk Download</h5>
          <Radio.Group
            options={[
              {
                label: 'wget',
                value: 'wget',
              },
              {
                label: 'curl',
                value: 'curl',
              }
            ]}
            value={bulkTool}
            onChange={(e) => {
              setBulkTool(e.target.value)
              setBulkText(BulkDownload(e.target.value, stid))
            }}
            optionType='button'
          />
        </Space>
        <TextArea
          placeholder='Bulk download commands'
          value={bulkText}
          autoSize={{
            minRows: 2,
            maxRows: 6,
          }}
        />
      </ConfigProvider>
    </>
  )
}
export default DatasetDetails