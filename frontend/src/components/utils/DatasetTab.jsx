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

const DatasetTab = () => {
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
        url: 'http://192.168.60.116:5001/datasets',
      })
        .then((response) => {
          setDataSrc(response.data.data)
          setDataCol(response.data.attributes)
          setFresh(false)
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
