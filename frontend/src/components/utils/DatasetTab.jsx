import React, { useEffect, useState } from 'react'
import { Tabs, Button, Space, Spin } from 'antd'
import DatasetList from './DatasetList'
import DatasetTable from './DatasetTable'
import DatasetStatic from './DatasetStatic'
import {
  UnorderedListOutlined,
  BarChartOutlined,
  TableOutlined,
  DownloadOutlined,
  ReloadOutlined,
  DotNetOutlined,
  ApartmentOutlined,
} from '@ant-design/icons'
import loadingTips from '../charts/LoadingTip'
import Loading from '../charts/Loading'
import { useQuery, QueryCache } from 'react-query'
import axios from 'axios'
import DatasetGraph from './DatasetGraph'

const {enterLoading, quitLoading} = Loading

const DatasetTab = (props) => {
  const { sendData } = props // send src metadata to parent
  const [dataSrc, setDataSrc] = useState([])
  const [dataCol, setDataCol] = useState({})
  const [fresh, setFresh] = useState(true)
  const [loadings, setLoadings] = useState([])
  const [currTip, setCurrTip] = useState(loadingTips[0])
  const items = [
    {
      key: '1',
      label: (
        <span>
          <UnorderedListOutlined />
          List
        </span>
      ),
      children: <DatasetList src={dataSrc} col={dataCol} />,
    },
    {
      key: '2',
      label: (
        <span>
          <TableOutlined />
          Table
        </span>
      ),
      children: <DatasetTable src={dataSrc} col={dataCol} />,
    },
    {
      key: '3',
      label: (
        <span>
          <ApartmentOutlined />
          Graph
        </span>
      ),
      children: <DatasetGraph config={{
        containerId: "viz",
        neo4j: {
            serverUrl: "neo4j://localhost:7687",
            serverUser: "neo4j",
            serverPassword: "biorzh123456",
        },
        labels: {
          ST: {
              caption: 'id',  
          },
          SC: {
              caption: 'id',   
          },
        },
        relationships: {
          PAIR:{
              thickness: 1,  //String：线段粗细，用作边缘厚度的属性名。默认为1。
              caption: true,  //Boolean：如果设置为true，则关系类型将显示为边缘标题。或String：用作边缘标题的属性名。
              font: { size: 12, color: '#606266' }  // 关系节点文字大小颜色
          },
        },
        initialCypher: "MATCH (n)-[r:INTERACTS]->(m) RETURN *"
      }} />,
    },
    {
      key: '4',
      label: (
        <span>
          <BarChartOutlined />
          Static
        </span>
      ),
      children: <DatasetStatic src={dataSrc} col={dataCol} />,
    },
  ]

  const response = useQuery({
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


  useEffect(() => {
    if(response.status === 'success' && typeof response.data !== 'undefined' && fresh) {
      setCurrTip(loadingTips[1])
      setDataSrc(response.data.data)
      setDataCol(response.data.attributes)
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
      setFresh(false)
    }
  }, [response.data, fresh])

  const operations = (
    <Space>
      <Button
        onClick={() => {
          enterLoading(0,setLoadings)
          setCurrTip(loadingTips[0])
          response.refetch().then(()=>{
            setFresh(true)
          })
        }}
        icon={<ReloadOutlined />}>
        Reload
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
    <Tabs defaultActiveKey="1" items={items} tabBarExtraContent={operations} />
    </Spin>
  )
}
export default DatasetTab
