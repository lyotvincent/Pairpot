import React, { useEffect, useState } from 'react'
import { Tabs, Button, Space } from 'antd'
import DatasetList from './DatasetList'
import DatasetTable from './DatasetTable'
import DatasetStatic from './DatasetStatic'
import {
  UnorderedListOutlined,
  BarChartOutlined,
  TableOutlined,
  DownloadOutlined,
  ReloadOutlined,
} from '@ant-design/icons'
import axios from 'axios'

const DatasetTab = (props) => {
  const { sendData } = props // send src metadata to parent
  const [dataSrc, setDataSrc] = useState([])
  const [dataCol, setDataCol] = useState({})
  const [fresh, setFresh] = useState(false)
  const items = [
    {
      key: '1',
      label: (
        <span>
          <UnorderedListOutlined />
          List
        </span>
      ),
      children: <DatasetList src={dataSrc} col={dataCol} />,
    },
    {
      key: '2',
      label: (
        <span>
          <TableOutlined />
          Table
        </span>
      ),
      children: <DatasetTable src={dataSrc} col={dataCol} />,
    },
    {
      key: '3',
      label: (
        <span>
          <BarChartOutlined />
          Static
        </span>
      ),
      children: <DatasetStatic src={dataSrc} col={dataCol} />,
    },
  ]

  const operations = (
    <Space>
      <Button
        onClick={() => {
          console.log('Reload')
          setFresh(true)
        }}
        icon={<ReloadOutlined />}>
        Reload
      </Button>
      <Button icon={<DownloadOutlined />}>Download</Button>
    </Space>
  )
  useEffect(() => {
    if (dataSrc.length === 0 || fresh) {
      console.log('loading data')
      axios({
        method: 'GET',
        url: 'http://localhost:5001/datasets',
      })
        .then((response) => {
          setDataSrc(response.data.data)
          setDataCol(response.data.attributes)
          setFresh(false)
          // calculate metadata and send to Database
          let src = response.data.data
          let attr = response.data.attributes
          let attrIndex = {}
          for (let i = 0; i < attr.length; i++) {
            let colname = attr[i]
            attrIndex[colname] = i
          }
          let filterName = ['species', 'tissues', 'technologies']
          let filterItem = filterName.map((fc) => {
            let Ori_col = Array.from(
              new Set(src.map((item) => item[attrIndex[fc]]))
            )
            let Uni_col = []
            for (let col of Ori_col) {
              if (col !== null) {
                Uni_col = [...Uni_col, ...col.split(';')]
              }
            }
            Uni_col = Array.from(new Set(Uni_col)).sort()
            Uni_col = Uni_col.map((item) => ({ text: item, value: item }))
            return Uni_col
          })
          let filterCol = {}
          for (let i = 0; i < filterName.length; i++) {
            filterCol[filterName[i]] = filterItem[i]
          }
          sendData(filterCol)
          //console.log(filterCol['technologies'].map((item) => item.text))
        })
        .catch((error) => {
          console.log(error)
        })
    }
  }, [fresh])
  return (
    <Tabs defaultActiveKey="1" items={items} tabBarExtraContent={operations} />
  )
}
export default DatasetTab
