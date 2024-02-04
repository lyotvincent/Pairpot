import React from 'react'
import { Tabs, Button } from 'antd'
import SampleList from './SampleList'
import SampleTable from './SampleTable'
import {
  UnorderedListOutlined,
  BarChartOutlined,
  TableOutlined,
  DownloadOutlined,
} from '@ant-design/icons'
const onChange = (key) => {
  console.log(key)
}
const operations = (
  <Button>
    <DownloadOutlined />
    Download
  </Button>
)
const items = [
  {
    key: '1',
    label: (
      <span>
        <UnorderedListOutlined />
        List
      </span>
    ),
    children: <SampleList />,
  },
  {
    key: '2',
    label: (
      <span>
        <TableOutlined />
        Table
      </span>
    ),
    children: <SampleTable />,
  },
  {
    key: '3',
    label: (
      <span>
        <BarChartOutlined />
        Static
      </span>
    ),
    children: 'Content of Tab Pane 3',
  },
]
const SampleTab = () => (
  <Tabs
    defaultActiveKey="1"
    items={items}
    onChange={onChange}
    tabBarExtraContent={operations}
  />
)
export default SampleTab
