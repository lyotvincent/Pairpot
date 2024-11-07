import React, { useState, useContext, useEffect } from 'react'
import '../App.css'
import {
  HomeOutlined,
  SearchOutlined,
  ToolOutlined,
  QuestionCircleOutlined,
  DatabaseOutlined,
  CloudUploadOutlined,
  CloudDownloadOutlined,
  PaperClipOutlined,
} from '@ant-design/icons'
import { Menu, ConfigProvider, Layout, Radio, theme, Space, Divider, Card } from 'antd'
import { NavLink } from 'react-router-dom'
import { useLocation } from 'react-router-dom'
import useToken from 'antd/es/theme/useToken'
import logoFig from "../assets/img/mylogo.png"
import Loading from './charts/Loading'
import { useNavigate } from 'react-router-dom'
import axios from 'axios'

const { Header } = Layout



const Nav = () => {
  const [themeToken, setThemeToken] = useState("secret")
  const [key, setKey] = useState("home")
  const location = useLocation()
  const { token } = theme.useToken()
  const handleThemeChange = (e) => {
    setThemeToken(e.target.value)
  }
  const navigate = useNavigate()

  const items = [
    {
      label: (
        <NavLink to="" style={{ textDecorationLine: 'none' }}>
          Home
        </NavLink>
      ),
      key: '/home',
      icon: <HomeOutlined />,
    },
    // {
    //   label: (
    //     <NavLink to="/browse" style={{ textDecorationLine: 'none' }}>
    //       Browse
    //     </NavLink>
    //   ),
    //   key: '/browse',
    //   icon: <SearchOutlined />,
    // },
    {
      label: (
        <NavLink to="/browse" style={{ textDecorationLine: 'none' }}>
          Browse
        </NavLink>
      ),
      key: '/browse',
      icon: <SearchOutlined />,
    },
    {
      label: (
        <NavLink to="/database" style={{ textDecorationLine: 'none' }}>
          Database
        </NavLink>
      ),
      key: '/database',
      icon: <DatabaseOutlined />,
    },
    // {
    //   label: (
    //     <NavLink to="/tools" style={{ textDecorationLine: 'none' }}>
    //       Tools
    //     </NavLink>
    //   ),
    //   key: 'tools',
    //   icon: <ToolOutlined />,
    // },
    // {
    //   label: 'Download',
    //   key: '/download',
    //   icon: <CloudDownloadOutlined />,
    // },
    // {
    //   label: (
    //     <NavLink to={`/submit`} style={{ textDecorationLine: 'none' }}>
    //       Submit
    //     </NavLink>
    //   ),
    //   key: 'submit',
    //   icon: <CloudUploadOutlined />,
    // },
    {
      label: (
        <NavLink to={`/help`} style={{ textDecorationLine: 'none' }}>
          Help
        </NavLink>
      ),
      key: '/help',
      icon: <QuestionCircleOutlined />,
    },
    // {
    //   label: (
    //     <NavLink to={`/news`} style={{ textDecorationLine: 'none' }}>
    //       News
    //     </NavLink>
    //   ),
    //   key: '/News',
    //   icon: <PaperClipOutlined />,
    // },
  ]

  useEffect(() => {
    if (location.pathname !== key) {
      setKey(location.pathname)
    }
  }, [location])

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


        <div className="demo-logo" style={{ width: 200, marginLeft: 5, marginRight: 5, backgroundImage: `url(${logoFig})`, backgroundSize: 'cover', backgroundPosition: 'center' }} />

        <Menu
          mode="horizontal"
          style={{
            flex: 1,
          }}
          items={items}
          selectedKeys={key}
          defaultSelectedKeys="home"
          onSelect={(e) => {
            setKey(e.key)
          }}
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
