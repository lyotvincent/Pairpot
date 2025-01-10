import React, { useEffect, useState, useImperativeHandle } from 'react'
import { Tabs, Button, Space, Spin } from 'antd'
import DatasetList from './DatasetList'
import DatasetTable from './DatasetTable'
import DatasetStatic from './DatasetStatic'
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
import { useQuery, QueryCache } from 'react-query'
import axios from 'axios'
import DatasetGraph from './DatasetGraph'

const {enterLoading, quitLoading} = Loading

const DatasetTab = (props) => {
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
  const OriginResponse = useQuery({
    queryKey: ['db'],
    queryFn: () => axios.get('/api/datasets').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })
  const items = [
    {
      key: '1',
      label: (
        <span>
          <UnorderedListOutlined />
          List
        </span>
      ),
      children: <DatasetList oriSrc={OriginResponse.data?.data} src={dataSrc} col={dataCol} />,
    },
    {
      key: '2',
      label: (
        <span>
          <TableOutlined />
          Table
        </span>
      ),
      children: <DatasetTable src={dataSrc} col={dataCol} oriSrc={OriginResponse.data?.data} />,
    },
    {
      key: '3',
      label: (
        <span>
          <ApartmentOutlined />
          Graph
        </span>
      ),
      children: <DatasetGraph config={graphCfg} src={dataSrc} col={dataCol}/>,
    },
    {
      key: '4',
      label: (
        <span>
          <BarChartOutlined />
          Static
        </span>
      ),
      children: <DatasetStatic 
      src={dataSrc} 
      col={dataCol} 
      _label={dataLabel} 
      _item={dataItem}
      onTab={tabKey === '4' ? true : false} />,
    },
  ]


  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Fresh": setFresh, // Trigger for useEffect
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
    "GraphConfig": setGraphCfg, // set GraphConfigs
  }))

  useEffect(() => {
    // console.log(response.data)
    // console.log(dataLabel)
    // console.log(dataItem)
    if(response.status === 'success' && typeof response.data !== 'undefined' && fresh) {
      setCurrTip(loadingTips[1])
      enterLoading(0, setLoadings)
      setDataSrc(response.data.data)
      setDataCol(response.data.attributes)
      setdataLabel(response.data.label) // add selected label
      setdataItem(response.data.item) // add selected item
      setFresh(false)
      // calculate metadata and send to Database
      let src = response.data.data
      let attr = response.data.attributes
      let attrIndex = {}
      for (let i = 0; i < attr.length; i++) {
        let colname = attr[i]
        attrIndex[colname] = i
      }
      let filterName = ['species', 'tissues', 'technologies']
      let filterItem = filterName.map((fc) => {
        let Ori_col = Array.from(
          new Set(src.map((item) => item[attrIndex[fc]]))
        )
        let Uni_col = []
        for (let col of Ori_col) {
          if (col !== null) {
            Uni_col = [...Uni_col, ...col.split(';')]
          }
        }
        Uni_col = Array.from(new Set(Uni_col)).sort()
        Uni_col = Uni_col.map((item) => ({ text: item, value: item }))
        return Uni_col
      })
      let filterCol = {}
      for (let i = 0; i < filterName.length; i++) {
        filterCol[filterName[i]] = filterItem[i]
      }
      sendData(filterCol)
      quitLoading(0, setLoadings)
    }
  }, [response.data, fresh])

  useEffect(() => {
    if(showAll && OriginResponse.status === 'success' && typeof OriginResponse.data !== 'undefined'){
      setCurrTip(loadingTips[1])
      setDataSrc(OriginResponse.data.data)
      setDataCol(OriginResponse.data.attributes)
      setdataLabel("all")
      setdataItem(OriginResponse.data.item)
      quitLoading(0, setLoadings)
      setShowAll(false)
    }
  }, [showAll, OriginResponse.data])

  const operations = (
    <Space>
      {/* <Button 
      onClick={() => {
        enterLoading(0,setLoadings)
        setCurrTip(loadingTips[0])
        setShowAll(true)
      }}
      icon={<ProductOutlined />}>
        Display All
      </Button> */}
      <Button
        onClick={() => {
          enterLoading(0,setLoadings)
          setCurrTip(loadingTips[0])
          response.refetch().then(()=>{
            setFresh(true)
          })
        }}
        icon={<ReloadOutlined />}>
        ReFresh
      </Button>
      {/* <Button icon={<DownloadOutlined />}>Download</Button> */}
    </Space>
  )

  useEffect(() => {
    setCurrTip(loadingTips[0])
  }, [])

  // useEffect(() => {
  //   enterLoading(0, setLoadings)
  //   setCurrTip(loadingTips[0])
  //   if (dataSrc.length === 0 || fresh) {
  //     axios({
  //       method: 'GET',
  //       url: '/api/datasets',
  //     })
  //       .then((response) => {
  //         setCurrTip(loadingTips[2])
  //         setDataSrc(response.data.data)
  //         setDataCol(response.data.attributes)
  //         setFresh(false)
  //         // calculate metadata and send to Database
  //         let src = response.data.data
  //         let attr = response.data.attributes
  //         let attrIndex = {}
  //         for (let i = 0; i < attr.length; i++) {
  //           let colname = attr[i]
  //           attrIndex[colname] = i
  //         }
  //         let filterName = ['species', 'tissues', 'technologies']
  //         let filterItem = filterName.map((fc) => {
  //           let Ori_col = Array.from(
  //             new Set(src.map((item) => item[attrIndex[fc]]))
  //           )
  //           let Uni_col = []
  //           for (let col of Ori_col) {
  //             if (col !== null) {
  //               Uni_col = [...Uni_col, ...col.split(';')]
  //             }
  //           }
  //           Uni_col = Array.from(new Set(Uni_col)).sort()
  //           Uni_col = Uni_col.map((item) => ({ text: item, value: item }))
  //           return Uni_col
  //         })
  //         let filterCol = {}
  //         for (let i = 0; i < filterName.length; i++) {
  //           filterCol[filterName[i]] = filterItem[i]
  //         }
  //         sendData(filterCol)
  //         quitLoading(0, setLoadings)
  //         //console.log(filterCol['technologies'].map((item) => item.text))
  //       })
  //       .catch((error) => {
  //         console.log(error)
  //       })
  //   }
  // }, [fresh])

  return (<Spin spinning={loadings[0]} size="large" tip={currTip}>
    <Tabs activeKey={tabKey} 
      items={items} 
      tabBarExtraContent={operations}
      onChange={(e) => setTabKey(e)} />
    </Spin>
  )
}
export default DatasetTab
