import { useState, useEffect, View } from 'react'
import { Container, Row, Col, Tab, Nav } from 'react-bootstrap'
import TrackVisibility from 'react-on-screen'
import { LoadingButton } from './utils/LoadingButton'
import { SwitchButton } from './utils/SwitchButton'
import 'bootstrap/dist/css/bootstrap.min.css'
import { ReactComponent as Female } from '../assets/img/female.svg'
import { ReactComponent as FemaleOrgans } from '../assets/img/female-organs.svg'
import SideMenu from './utils/SideMenu'
import { Breadcrumb, Layout, Menu, theme, Space, Button, Card } from 'antd'
import {
  MenuFoldOutlined,
  MenuUnfoldOutlined,
  DatabaseOutlined,
  ProfileOutlined,
  OneToOneOutlined,
} from '@ant-design/icons'
import DatasetTab from './utils/DatasetTab'
import SampleTab from './utils/SampleTab'
const { Header, Content, Footer, Sider } = Layout

export const Database = () => {
  const [collapsed, setCollapsed] = useState(false)
  const [srcMeta, setSrcMeta] = useState({}) // get the meta data from childrens
  const [selectedKey, setSelectedKey] = useState('datasets')
  const [visible, setVisible] = useState(true)
  const onClickTopMenu = (e) => {
    setSelectedKey(e.key)
  }

  useEffect(() => {
    setVisible(true) // 默认显示第一个Content组件
  }, [])
  var topMenuItems = [
    {
      label: 'Datasets',
      key: 'datasets',
      icon: <DatabaseOutlined />,
    },
    {
      label: 'Samples',
      key: 'samples',
      icon: <ProfileOutlined />,
    },
    {
      label: 'Pairs',
      key: 'pairs',
      icon: <OneToOneOutlined />,
    },
  ]
  const {
    token: { colorBgContainer },
  } = theme.useToken()
  return (
    <Layout
      style={{
        minHeight: '100vh',
      }}>
      <Sider
        width={200}
        breakpoint="lg"
        trigger={null}
        collapsible
        collapsed={collapsed}>
        <SideMenu items={srcMeta} />
      </Sider>
      <Layout>
        <Header
          style={{
            padding: 0,
          }}>
          <Space>
            <Button
              type="text"
              icon={collapsed ? <MenuUnfoldOutlined /> : <MenuFoldOutlined />}
              onClick={() => setCollapsed(!collapsed)}
              style={{
                fontSize: '16px',
                width: 64,
                height: 64,
                color: '#dddddd',
              }}
            />
            <Menu
              width={'70%'}
              onClick={onClickTopMenu}
              selectedKeys={[selectedKey]}
              mode="horizontal"
              items={topMenuItems}
              theme="dark"
            />
          </Space>
        </Header>
        <Content
          style={{
            padding: '0rem 1rem',
          }}>
          {visible && selectedKey === 'datasets' && (
            <DatasetTab sendData={setSrcMeta} />
          )}
          {visible && selectedKey === 'samples' && <SampleTab />}
          {visible && selectedKey === 'pairs' && <div>pairs</div>}
        </Content>
        <Footer
          style={{
            textAlign: 'center',
          }}>
          pairST ©2023 Created by zhihan Ruan
        </Footer>
      </Layout>
    </Layout>
  )
}
