import { useEffect, useState, View } from 'react'
import { Container, Row, Col, ListGroup, ListGroupItem } from 'react-bootstrap'
import { Menu } from 'antd'
import {
  ExperimentOutlined,
  SwitcherOutlined,
  YuqueOutlined,
} from '@ant-design/icons'
import PropTypes from 'prop-types'

import '../theme/dark'
import '../theme/vintage'

function getItem(label, key, icon, children, type) {
  return {
    key,
    icon,
    children,
    label,
    type,
  }
}
const InitialItems = [
  getItem('Species', 'sub1', <YuqueOutlined />, [
    getItem('Homo sapiens(Human)', '1'),
    getItem('Mus musculus(Mouse)', '2'),
    getItem('Danio rerio(Zebrafish)', '3'),
    getItem('Macaca fascicularis (Monkeys)', '4'),
  ]),
  getItem('Tissues', 'sub2', <SwitcherOutlined />, [
    getItem('Brain', '5'),
    getItem('Embryo', '6'),
    getItem('Skin', '7'),
    getItem('Kidney', '8'),
  ]),
  getItem('Technologies', 'sub3', <ExperimentOutlined />, [
    getItem('10X Visium', '9'),
    getItem('ST', '10'),
    getItem('GeoMx DSP', '11'),
    getItem('Stereo-seq', '12'),
  ]),
]

const SideMenu = ({ theme, title, items }) => {
  const [menuItems, setMenuItems] = useState(InitialItems)
  const [openKeys, setOpenKeys] = useState(['sub1'])
  var rootSubmenuKeys = ['sub1', 'sub2', 'sub4']
  const onOpenChange = (keys) => {
    const latestOpenKey = keys.find((key) => openKeys.indexOf(key) === -1)
    if (latestOpenKey && rootSubmenuKeys.indexOf(latestOpenKey) === -1) {
      setOpenKeys(keys)
    } else {
      setOpenKeys(latestOpenKey ? [latestOpenKey] : [])
    }
  }
  const ItemIdx = {
    species: 0,
    tissues: 1,
    technologies: 2,
  }

  useEffect(() => {
    let __menuItems = InitialItems
    let usedKey = 0
    for (let key in items) {
      let idx = ItemIdx[key]
      __menuItems[idx].children = items[key].map((item, i) =>
        getItem(item.value, i + usedKey)
      )
      usedKey = usedKey + items[key].length
    }
    console.log(__menuItems)
    setMenuItems(__menuItems)
  }, [items])

  return (
    <Menu
      mode="inline"
      theme={theme}
      title={title}
      openKeys={openKeys}
      onOpenChange={onOpenChange}
      items={menuItems}
    />
  )
}

SideMenu.defaultProps = {
  theme: 'dark',
  title: 'Datasets',
  items: { Species: [], Tissues: [], Technologies: [] },
}

SideMenu.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  items: PropTypes.array,
}

export default SideMenu
