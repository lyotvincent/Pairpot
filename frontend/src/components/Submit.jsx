import { FileTextOutlined, LinkOutlined } from '@ant-design/icons'
import React, { useState, useEffect } from 'react'
import {
  Menu,
  Button,
  Input,
  Col,
  Form,
  InputNumber,
  Radio,
  Row,
  Select,
  Space,
  Layout,
  theme,
  Result,
  Modal,
} from 'antd'
import { useLocation, useNavigate } from 'react-router-dom'
import axios from 'axios'
import SubmitMeta from './utils/SubmitMeta'
import SubmitLink from './utils/SubmitLink'
const { Content, Sider } = Layout
const { useToken } = theme

const menuItems = [
  {
    key: 'link',
    label: 'Link Paired Data',
    icon: React.createElement(LinkOutlined),
  },
  {
    key: 'new',
    label: 'New SRT Data',
    icon: React.createElement(FileTextOutlined),
    children: [
      {
        key: 'new-meta',
        label: 'New Meta Data',
      },
      {
        key: 'new-expr',
        label: 'New expression file',
      },
      {
        key: 'new-correct',
        label: 'Correct Meta Data',
      },
    ],
  },
]

const Submit = () => {
  const navigate = useNavigate()
  const location = useLocation()
  const [submitForm] = Form.useForm()
  const [selectedKey, setSelectedKey] = useState('link')
  const [visible, setVisible] = useState(true)
  const onClickSideMenu = (e) => {
    setSelectedKey(e.key)
  }
  const [result, setResult] = useState({ show: false, state: '', message: '' })
  const onFinish = (values) => {
    console.log('Received values of form: ', values)
    axios
      .post('http://localhost:5522/submit', {
        data: values,
        headers: {
          'Content-Type': 'application/json',
        },
      })
      .then((response) => {
        setResult({
          show: true,
          state: response.data.state,
          message: response.data.message,
        })
      })
  }
  const { token } = useToken();
  useEffect(() => {
    setVisible(true) // 默认显示第一个Content组件
  }, [])

  return (
    <Layout>
      <Sider
        width={200}
        breakpoint="lg"
        trigger={null}
        collapsible
        style={{ position: 'sticky' }}>
        {' '}
        <Menu
          theme="dark"
          mode="inline"
          onClick={onClickSideMenu}
          defaultSelectedKeys={'sub1'}
          style={{
            height: '100%',
            borderRight: 0,
          }}
          items={menuItems}
        />
      </Sider>
      <Content style={{ padding: '20px 20px' }}>
        <div
          style={{
            padding: 24,
            minHeight: 380,
            background:  token.colorBgContainer,
            borderRadius: token.borderRadiusLG,
          }}>
          {visible && selectedKey === 'link' && (
            <SubmitLink
              onFinish={onFinish}
              location={location}
              navigate={navigate}
              submitForm={submitForm}
            />
          )}
          {visible && selectedKey === 'new-meta' && (
            <SubmitMeta onFinish={onFinish} />
          )}
        </div>
      </Content>
      <Modal
        open={result.show}
        onOk={() => {
          setResult({ show: false, state: '', message: '' })
        }}
        onCancel={() => {
          setResult({ show: false, state: '', message: '' })
        }}
        style={{ top: '15%' }}
        width={1000}>
        <Result
          status={result.state}
          title={result.message}
          subTitle="Cloud server configuration takes 1-5 minutes, please wait."
          extra={[
            <Button
              type="primary"
              key="gohome"
              onClick={() => {
                setResult({ show: false, state: '', message: '' })
                navigate('/')
              }}>
              Go HomePage
            </Button>,
            <Button
              key="submit"
              onClick={() => {
                setResult({ show: false, state: '', message: '' })
                navigate('/database')
              }}>
              Submit Again
            </Button>,
          ]}
        />
      </Modal>
    </Layout>
  )
}

export default Submit
