import React, { useState, useContext } from 'react'
import '../App.css'
import {
  HomeOutlined,
  SearchOutlined,
  ToolOutlined,
  QuestionCircleOutlined,
  DatabaseOutlined,
  CloudUploadOutlined,
  CloudDownloadOutlined,
} from '@ant-design/icons'
import { Menu, ConfigProvider, Layout, Radio, theme, Space, Divider, Card } from 'antd'
import { NavLink } from 'react-router-dom'
import useToken from 'antd/es/theme/useToken'
const { Header } = Layout

const items = [
  {
    label: (
      <NavLink to="" style={{ textDecorationLine: 'none' }}>
        Home
      </NavLink>
    ),
    key: 'home',
    icon: <HomeOutlined />,
  },
  {
    label: (
      <NavLink to="/browse" style={{ textDecorationLine: 'none' }}>
        Browse
      </NavLink>
    ),
    key: 'browse',
    icon: <SearchOutlined />,
  },
  {
    label: (
      <NavLink to="/database" style={{ textDecorationLine: 'none' }}>
        Database
      </NavLink>
    ),
    key: 'db',
    icon: <DatabaseOutlined />,
  },
  {
    label: (
      <NavLink to="/tools" style={{ textDecorationLine: 'none' }}>
        Tools
      </NavLink>
    ),
    key: 'tools',
    icon: <ToolOutlined />,
    // children: [
    //   {
    //     label: 'Lasso-Refine',
    //     key: 'lr',
    //   },
    //   {
    //     label: 'Spatial Variable Genes',
    //     key: 'svg',
    //   },
    //   {
    //     label: 'Lasso-Deconvolution',
    //     key: 'ld',
    //   },
    //   {
    //     label: 'Cell Communication',
    //     key: 'cc',
    //   },
    //   {
    //     label: 'Function Analysis',
    //     key: 'fa',
    //   },
    // ],
  },
  {
    label: 'Download',
    key: 'dl',
    icon: <CloudDownloadOutlined />,
  },
  {
    label: (
      <NavLink to={`/submit`} style={{ textDecorationLine: 'none' }}>
        Submit
      </NavLink>
    ),
    key: 'submit',
    icon: <CloudUploadOutlined />,
  },
  {
    label: (
      <NavLink to={`/contact`} style={{ textDecorationLine: 'none' }}>
        Help
      </NavLink>
    ),
    key: 'help',
    icon: <QuestionCircleOutlined />,
  },
]

const Nav = () => {
  const [themeToken, setThemeToken] = useState("secret")
  const { token } = theme.useToken()
  const handleThemeChange = (e) => {
    setThemeToken(e.target.value)
  }
  return (
    <ConfigProvider
      theme={{
        components: {
          Menu: {
            iconSize: 24,
            fontSize: 24,
            colorBgContainer: '#ddf7d4',
            iconMarginInlineEnd: 2,
            colorInfoTextHover: token.colorPrimary,
            itemHoverColor: token.colorPrimary,
          }
        },
      }}>
      <Header
        style={{
          display: 'flex',
          alignItems: 'center',
          backgroundColor: '#ddf7d4',
          height: 60,
          position: 'sticky',
          top: 0,
          zIndex: 1060,
          width: '100%',
        }}>


        <div className="demo-logo" style={{ width: 150, backgroundColor: '#cdf3bf', margin: 0 }} />

        <Menu
          mode="horizontal"
          style={{
            flex: 1,
          }}
          items={items}
        />
        {/* <Space>
          <span>Theme: </span>
          <Radio.Group
            value={themeToken}
            onChange={handleThemeChange}
            defaultValue="bright"
            buttonStyle="solid"
            size="small"
            style={{ marginRight: 20 }}>
            <Radio.Button value="bright">Bright</Radio.Button>
            <Radio.Button value="secret">Secret</Radio.Button>

          </Radio.Group>
        </Space> */}
      </Header>
    </ConfigProvider>
  )
}
export default Nav
