import { useState, View } from 'react'
import { Menu } from 'antd'
import {
  DatabaseOutlined,
  OneToOneOutlined,
  OrderedListOutlined,
  ProfileOutlined,
  TableOutlined,
} from '@ant-design/icons'
import PropTypes from 'prop-types'

const items = [
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
const TopMenu = () => {
  const [current, setCurrent] = useState('statistic')
  const onClick = (e) => {
    console.log('click ', e)
    setCurrent(e.key)
  }
  return (
    <Menu
      width={'70%'}
      onClick={onClick}
      selectedKeys={[current]}
      mode="horizontal"
      items={items}
      theme="dark"
    />
  )
}
export default TopMenu
