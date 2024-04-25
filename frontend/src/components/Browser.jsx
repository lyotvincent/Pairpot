import { DotChartOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense } from 'react'
import {
  Menu,
  Button,
  Input,
  Col,
  Form,
  Steps,
  Layout,
  theme,
  Result,
  Modal,
  Divider
} from 'antd'
import { useLocation, useNavigate } from 'react-router-dom'
import axios from 'axios'
const { Content, Sider } = Layout

const sideMenuItems = [
  {
    label: 'MetaInfo',
    key: 'MetaInfo',
    icon: <FileSearchOutlined />
  },
  {
    label: 'Layer View',
    key: 'Layer View',
    icon: <DotChartOutlined />
  },
  {
    label: 'Pair View',
    key: 'Pair View',
    icon: <OneToOneOutlined />
  },
  {
    label: 'Spatial Markers',
    key: 'Spatial Markers',
    icon: <HeatMapOutlined />,
  },
  {
    label: 'GO Enrichment',
    key: 'GO Enrichment',
    icon: <HeatMapOutlined />,
  },
  {
    label: 'KEGG Enrichment',
    key: 'KEGG Enrichment',
    icon: <HeatMapOutlined />,
  },
  {
    label: 'Cell Interactions',
    key: 'Cell Interactions',
    icon: <HeatMapOutlined />,
  },
]

const MetaInfo = lazy(() => (import('./utils/DatasetDetails')))
const LayerView = lazy(() => (import('./charts/SpScatter')))

const Browser = () => {
  const navigate = useNavigate()
  const location = useLocation()
  const [collapsed, setCollapsed] = useState(false)
  const [selectedKey, setSelectedKey] = useState(sideMenuItems[0].key)
  const onClickSideMenu = (e) => {
    setSelectedKey(e.key)
  }
  return (
    <Layout>
      <Sider
        width={200}
        breakpoint="lg"
        collapsible collapsed={collapsed} onCollapse={(value) => setCollapsed(value)}>
        <Menu
          theme="dark"
          mode="inline"
          onClick={onClickSideMenu}
          defaultSelectedKeys={'sub1'}
          selectedKeys={selectedKey}
          style={{
            width: collapsed ? 78 : 198,
            height: '100%',
            borderRight: 0,
            position: 'fixed'
          }}
          items={sideMenuItems}
        />
      </Sider>
      <Content style={{ padding: '20px 20px' }}>
        <Suspense fallback={<h1>Loading...</h1>}>
          <MetaInfo
            location={location}
            navigate={navigate} />
        </Suspense>
        <Divider />
        <Suspense fallback={<h1>Loading...</h1>}>
          <LayerView />
        </Suspense>
        {/* <div style={{color: 'black'}}>{JSON.stringify(location.state)}</div> */}
      </Content>
    </Layout>
  )
}
export default Browser