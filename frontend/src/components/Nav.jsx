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
import { Menu, ConfigProvider, Layout } from 'antd'
import { NavLink } from 'react-router-dom'
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
    label: 'Tools',
    key: 'tools',
    icon: <ToolOutlined />,
    children: [
      {
        label: 'Lasso-Refine',
        key: 'lr',
      },
      {
        label: 'Spatial Variable Genes',
        key: 'svg',
      },
      {
        label: 'Lasso-Deconvolution',
        key: 'ld',
      },
      {
        label: 'Cell Communication',
        key: 'cc',
      },
      {
        label: 'Function Analysis',
        key: 'fa',
      },
    ],
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
  {
    label: 'About',
    key: 'about',
  },
]

const Nav = () => {
  return (
    <Header
      style={{
        display: 'flex',
        alignItems: 'center',
        backgroundColor: '#100c2a',
        height: 70,
        position: 'sticky',
        top: 0,
        zIndex: 1060,
        width: '100%',
        display: 'flex',
      }}>
      <div className="demo-logo" />
      <ConfigProvider
        theme={{
          components: {
            Menu: {
              iconSize: 20,
              fontSize: 20,
              darkItemSelectedBg: '#000000',
              horizontalItemSelectedColor: '#001529',
              colorBgBase: '#100c2a',
              iconMarginInlineEnd: 8,
            },
          },
        }}>
        <Menu
          mode="horizontal"
          theme="dark"
          style={{
            flex: 1,
            minWidth: 0,
            backgroundColor: '#100c2a',
            color: 'white',
          }}
          items={items}
          selectable={0}
        />
      </ConfigProvider>
    </Header>
  )
}
export default Nav
