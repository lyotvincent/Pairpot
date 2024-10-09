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

import Anndata2sce from '../utils/ToolsCode/Anndata2sce'
import Anndata2seurat from '../utils/ToolsCode/Anndata2seurat'


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
    {
      key: '1',
      label: (
        <span>
          Anndata to SingleCellExperiment
        </span>
      ),
      children: <Anndata2sce />
    },
    {
      key: '2',
      label: (
        <span>
          Anndata to Seurat
        </span>
      ),
      children: <Anndata2seurat />
    },
  ]


  return (
    <Tabs activeKey={tabKey}
      items={items}
      // tabBarExtraContent={operations}
      onChange={(e) => setTabKey(e)} />
  )
}
export default HelpTab
