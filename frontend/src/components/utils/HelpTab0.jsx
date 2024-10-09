import React, { useEffect, useState, useImperativeHandle, Children } from 'react'
import { Tabs, Button, Space, Spin } from 'antd'
import {
  UnorderedListOutlined,
  BarChartOutlined,
  TableOutlined,
  // DownloadOutlined,
  ReloadOutlined,
  // DotNetOutlined,
  ProductOutlined,
  ApartmentOutlined,
} from '@ant-design/icons'
import loadingTips from '../charts/LoadingTip'
import Loading from '../charts/Loading'

import CellChat from "./ToolsCode/CellChat"
import CellTrek from './ToolsCode/CellTrek'
import ITALK from './ToolsCode/ITalk'
import Card from './ToolsCode/Card'
import SeuratDeconv from './ToolsCode/SeuratDeconv'
import SpaTalkDeconv from './ToolsCode/SpaTalkDeconv'
import Spacexr from './ToolsCode/Spacexr'

const { enterLoading, quitLoading } = Loading

const HelpTab = (props) => {
  const { sendData, response, onRef } = props // send src metadata to parent
  const [dataSrc, setDataSrc] = useState([])
  const [dataCol, setDataCol] = useState({})
  const [dataLabel, setdataLabel] = useState({})
  const [dataItem, setdataItem] = useState({})
  const [fresh, setFresh] = useState(true)
  const [tabKey, setTabKey] = useState("1")
  const [loadings, setLoadings] = useState([])
  const [showAll, setShowAll] = useState(false)
  const [currTip, setCurrTip] = useState(loadingTips[0])
  const [graphCfg, setGraphCfg] = useState(null)
  const items = [
    // {
    //   key: '1',
    //   label: (
    //     <span>
    //       <UnorderedListOutlined />
    //       CellChat
    //     </span>
    //   ),
    //   children: <CellChat/>
    // },
    // {
    //   key: '2',
    //   label: (
    //     <span>
    //       <TableOutlined />
    //       CellTrek
    //     </span>
    //   ),
    //   children: <CellTrek/>
    // },
    // {
    //   key: '3',
    //   label: (
    //     <span>
    //       <ApartmentOutlined />
    //       iTALK
    //     </span>
    //   ),
    //   children: <ITALK/>
    // },
    {
      key: '1',
      label: (
        <span>
          CARD
        </span>
      ),
      children: <Card />
    },
    {
      key: '2',
      label: (
        <span>
          RCTD
        </span>
      ),
      children: <Spacexr />
    },
    {
      key: '3',
      label: (
        <span>
          SpaTalk
        </span>
      ),
      children: <SpaTalkDeconv />
    },
    {
      key: '4',
      label: (
        <span>
          Seurat
        </span>
      ),
      children: <SeuratDeconv />
    }
  ]


  return (
    <Tabs activeKey={tabKey}
      items={items}
      // tabBarExtraContent={operations}
      onChange={(e) => setTabKey(e)} />
  )
}
export default HelpTab
