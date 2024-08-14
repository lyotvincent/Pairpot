import { DotChartOutlined, FileDoneOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, SelectOutlined, SlidersOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense, useRef,memo } from 'react'
import logoFig from "../assets/img/mylogo.png"
import {
  Menu,
  Button,
  Input,
  Col,
  Form,
  Steps,
  Layout,
  theme,
  Divider,
  Card,
  Switch,
  ConfigProvider,
  Anchor,
  Flex,
} from 'antd'
import ToggleAccordion from './utils/ToggleAccordion'
import { Link, useLocation, useNavigate } from 'react-router-dom'
import MarkerTable from './charts/MarkerTable'
import axios from 'axios'
import NetworkRelation from './charts/NetworkRelation'
import RelHeat from './charts/Relation-Heatmap-test'
import Search from './utils/Search'
import loadingTips from './charts/LoadingTip'
import Loading from './charts/Loading'
import { useMutation, useQuery, QueryCache } from "react-query"
const { Content, Sider, Header } = Layout
const {enterLoading} = Loading

const MetaInfo = memo(lazy(() => (import('./utils/DatasetDetails'))))
const LassoView = lazy(() => (import('./charts/ScScatter')))
const LayerView = lazy(() => (import('./charts/SpScatter')))
const PairView = lazy(() => (import('./charts/PairView')))

const Browser = () => {
  const navigate = useNavigate()
  const location = useLocation()
  const [Init, setInit] = useState(false)
  const [collapsed, setCollapsed] = useState(false)
  const [locState, setLocState] = useState(false)
  const [action, setAction] = useState(0)
  const commandRef = useRef()
  const onClickSideMenu = (e) => {

  }
  const fetchSPData = (state) => (
    axios.get('/api/query/sp', {
      params: {
        id: state.st?.dataset_id
      },
      responseType: 'blob',
    }).then((response) => {
      let spblob = response.data
      let _spfile = new File([spblob], "sp_meta.h5ad")
      return _spfile
    })
  )

  const fetchSCData = (state) => (
    axios.get('/api/query/sc', {
      params: {
        id: state.st?.dataset_id
      },
      responseType: 'blob',
    }).then((response) => {
      let scblob = response.data
      let _scfile = new File([scblob], "sc_sampled.h5ad")
      return _scfile
    })
  )

  const statusSP = useMutation({
    mutationKey: ['sp'],
    mutationFn: (state) => fetchSPData(state),
    staleTime: Infinity, 
    retry: false, 
    refetchOnWindowFocus: false, 
  })

  const statusSC = useMutation({
    mutationKey: ['sc'],
    mutationFn: (state) => fetchSCData(state),
    staleTime: Infinity, 
    retry: false, 
    refetchOnWindowFocus: false, 
  })

  const example = useQuery({
    queryKey: ['example'],
    queryFn: ()=> axios.get('/api/example',{
        params: {
          id: "STDS0000235"
        },
      }).then((response) => {
        let dataCol = response.data.attributes
        let spitem = response.data.data[0]
        let values = Object.fromEntries(
          dataCol.map((k, i) => [k, spitem[i]])
        )
        let scitem = response.data.data[1]
        let state = {
          st: values
        }
        if (typeof scitem !== 'undefined') {
          let scvalues = Object.fromEntries(
            dataCol.map((k, i) => [k, scitem[i]])
          )
          state['sc'] = scvalues
        }
        return state
      })
    ,
    staleTime: Infinity, 
    retry: false, 
    refetchOnWindowFocus: false, 
  })

  const SearchAnchor = React.createRef() // location of Search
  const MetaAnchor = React.createRef() // location of MetaAnchor
  const LassoAnchor = React.createRef() // location of LassoView
  const LayerAnchor = React.createRef() // location of LassoView
  const PairAnchor = React.createRef() // location of LassoView
  const MarkerAnchor = React.createRef() // location of LassoView
  const CPDBAnchor = React.createRef() // location of LassoView
  const sideMenuItems = [
    {
      label: 'Search',
      ref: SearchAnchor,
      key: 'Search',
      icon: <FileSearchOutlined />
    },
    {
      label: 'MetaInfo',
      ref: MetaAnchor,
      key: 'MetaInfo',
      icon: <FileDoneOutlined />
    },
    {
      label: 'LassoView',
      ref: LassoAnchor,
      key: 'LassoView',
      icon: <SelectOutlined />
    },
    {
      label: 'Layer View',
      ref: LayerAnchor,
      key: 'LayerView',
      icon: <DotChartOutlined />
    },
    {
      label: 'Pair View',
      ref: PairAnchor,
      key: 'PairView',
      icon: <OneToOneOutlined />
    },
    {
      label: 'Spatial Markers',
      ref: MarkerAnchor,
      key: 'SpatialMarkers',
      icon: <SlidersOutlined />,
    },
    {
      label: 'Cell Interactions',
      ref: CPDBAnchor,
      key: 'CellInteractions',
      icon: <HeatMapOutlined />,
    },
  ]
  const [selectedKey, setSelectedKey] = useState(sideMenuItems[0].key)
  const MetaRef = React.createRef()   // get the trigger function from MetaInfo
  const LassoRef = React.createRef()  // get the trigger function from LassoView
  const LayerRef = React.createRef()  // get the trigger function from LayerView
  const PairRef = React.createRef()   // get the trigger function from PairView
  const MarkerRef = React.createRef()   // get the trigger function from Marker Table
  const NetRef = React.createRef()   // get the trigger function from Network Table
  const CPDBRef = React.createRef()   // get the trigger function from CPDB heatmap

  const [componentLoad, setComponentLoad] = useState({  // whether component is called to load
    "LassoView": false,
    "LayerView": false,
    "PairView": false,
    "MarkerTable": false,
    "NetWorkRelation": false,
    "CPDBHeatmap": false,
  })

  useEffect(()=>{
    if(example.status === 'success' && typeof example.data !== 'undefined'){
      MetaRef.current?.Trigger(example.data)
      setLocState(example.data)
    }
  }, [example.data, MetaRef])

  useEffect(() => {
    if(typeof locState !== 'undefined'){
      statusSC.mutate(locState)
      statusSP.mutate(locState)
    }
  }, [locState])

  useEffect(() => {
    if (Init && location.state !== null) {
      setComponentLoad({  // reset ComponentLoad
        "LassoView": false,
        "LayerView": false,
        "PairView": false,
        "MarkerTable": false,
        "NetWorkRelation": false,
        "CPDBHeatmap": false,
      })

      // set loading status for each components
      LassoRef.current.Tip(loadingTips[3])
      enterLoading(0, LassoRef.current.Loading)
      LayerRef.current.Tip(loadingTips[3])
      enterLoading(0, LayerRef.current.Loading)
      PairRef.current.Tip(loadingTips[3])
      enterLoading(1, PairRef.current.Loading)
      MarkerRef.current.Tip(loadingTips[3])
      enterLoading(0, MarkerRef.current.Loading)
      enterLoading(1, MarkerRef.current.Loading)
      CPDBRef.current.Tip(loadingTips[3])
      enterLoading(0, CPDBRef.current.Loading)
      NetRef.current.Tip(loadingTips[3])
      enterLoading(0, NetRef.current.Loading)

      // mutation of datasets
      statusSC.mutate(location.state)
      statusSP.mutate(location.state)
    }
    else {
        setInit(true)
    }
  }, [location])

  useEffect(() => {
    if (statusSP.status === 'success' &&
      typeof statusSP.data !== "undefined") {
      if(!componentLoad['LayerView']){
        LayerRef.current?.Trigger("Reload")  // only reload once
      }
      if(!componentLoad['MarkerTable']){
        MarkerRef.current?.Trigger("Reload")
      }  
    }
    if (statusSP.status === 'success' &&
      statusSC.status === 'success' &&
      typeof statusSC.data !== "undefined" &&
      typeof statusSP.data !== "undefined") {
      if(!componentLoad['LassoView']){
        LassoRef.current?.Trigger("Reload")
      }
      if(!componentLoad['PairView']) {
        PairRef.current?.Trigger("Reload")
      }
      if(!componentLoad['NetWorkRelation']) {
        NetRef.current?.Trigger("Reload")
      }
      if(!componentLoad['CPDBHeatmap']){
        CPDBRef.current?.Trigger("Reload")
      } 
    }
  }, [statusSP.variables, statusSC.variables, PairRef])

  return (
    <ConfigProvider
      theme={{
        components: {
          Anchor: {
            linkPaddingInlineStart: 0
          },
          Menu: {
            iconSize: 18,
          }
        },
      }}
    >
      <Layout>
        <Sider
          width={200}
          breakpoint="lg"
          collapsible collapsed={collapsed} onCollapse={(value) => setCollapsed(value)}>
          <Menu
            mode="inline"
            onClick={onClickSideMenu}
            defaultSelectedKeys={'sub1'}
            selectedKeys={selectedKey}
            style={{
              width: collapsed ? 60 : 200,
              height: '100%',
              borderRight: 0,
              position: 'fixed'
            }}
          >
            {sideMenuItems.map((item) => (
              <Menu.Item
                key={item.key}
                icon={item.icon}
                onClick={(e) => {
                  setSelectedKey(e.key)
                  item.ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                }}>
                {item.label}
              </Menu.Item>
            ))}
          </Menu>
        </Sider>
        <Layout>
          <Content style={{ padding: '20px 20px' }}>
            <Col span={6} offset={8} id='Search'>
              <div ref={SearchAnchor} style={{
                width: '300px',
                height: "70px",
                marginBottom: 30,
                backgroundImage: `url(${logoFig})`,
                objectFit: "cover",
                backgroundSize: 'cover',
                alignItems: "center"
              }} />
            </Col>
            <Col span={18} offset={3}> <Search /> </Col>
            <Divider />
            {/* <div>compLoad:{JSON.stringify(componentLoad)}</div> */}
            <Suspense fallback={<h1>Loading for MetaInfo...</h1>}>
              <Card id="MetaInfo" ref={MetaAnchor}>
                <MetaInfo onRef={MetaRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LassoView...</h1>}>
              <Card>
                <ToggleAccordion header={
                  <h3 id="LassoView" ref={LassoAnchor}>{"Lasso View "}
                  </h3>}>
                  {<div>
                    (left) These are the assignments of each spot to clusters by an
                    automated clustering algorithm in a single slice. The clustering groups together spots that
                    have similar expression profiles. In this plot, spots are colored according
                    to their cluster assignment and projected on to the tissue image. Only spots
                    under tissue are used in the clustering algorithm.
                    <br />
                    (right) These are the assignments of each spot to clusters by an
                    automated clustering algorithm in a single slice. The clustering groups together spots that
                    have similar expression profiles. In this plot, spots are colored according
                    to their cluster assignment and projected on to the tissue image. Only spots
                    under tissue are used in the clustering algorithm.
                  </div>}
                </ToggleAccordion>
                <Button type="primary"
                  size='middle'
                  onClick={() => {
                    LassoAnchor.current.scrollIntoView({ behavior: 'smooth', block: 'start' })
                    setTimeout(() => { LassoRef.current.Tour(true) }, 600)
                  }}>
                  {"Begin Tour >>>"}
                </Button>
                <LassoView spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} location={location} onRef={LassoRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LayerView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="LayerView" ref={LayerAnchor}>Layer View</h3>}>
                  {<div>
                    (left) These are the assignments of each cell to clusters by an
                    automated clustering algorithm in single-cell or spatal transcriptomics data. The clustering groups together cells that
                    have similar expression profiles. These clusters are further annotated manually according to their marker genes for single-cell data.
                    In this plot, cells are colored according to their cluster assignment and projected to UMAP space. In this space, pairs of cells that are close to each other have
                    more similar gene expression profiles than cells that are distant from each other.
                    <br />
                    (right)In this plot, cells are not colored but projected to UMAP space.
                    After lasso selection in the left figure, the selected region series would display in this plot.
                    After Online refine, the refined region series would also display in this plot.
                  </div>}
                </ToggleAccordion>
                <Button
                  type="primary"
                  size='middle'
                  onClick={() => {
                    LayerAnchor.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                    setTimeout(() => {
                      LayerRef.current.Tour(true)
                    }, 600)
                  }}>
                  {"Begin Tour >>>"}
                </Button>
                <LayerView id='LayerView' query={true} spfile={statusSP.data} setCompLoad={setComponentLoad} onRef={LayerRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for PairView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="PairView" ref={PairAnchor}>Pair View</h3>}>
                  {<div>
                    (left) These are the assignments of each cell to clusters by an
                    automated clustering algorithm in single-cell data. The clustering groups together cells that
                    have similar expression profiles. These clusters are further annotated manually according to their marker genes.
                    In this plot, cells are downsampled to 3000 and colored according to their cluster assignment and projected to UMAP space. In this space, pairs of cells that are close to each other have
                    more similar gene expression profiles than cells that are distant from each other.
                    <br />
                    (right) Spots are colored by clustering assignment and shown in t-SNE space.
                    The axes correspond to the 2-dimensional embedding produced by the t-SNE
                    algorithm. In this space, pairs of spots that are close to each other have
                    more similar gene expression profiles than spots that are distant from each other.
                  </div>}
                </ToggleAccordion>
                <Button
                  type="primary"
                  size='middle'
                  onClick={() => {
                    PairAnchor.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                    setTimeout(() => {
                      PairRef.current.Tour(true)
                    }, 600)
                  }}>
                  {"Begin Tour >>>"}
                </Button>
                <PairView spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} location={location} onRef={PairRef} />
              </Card>
            </Suspense>
            <Divider />
            <Card>
              <ToggleAccordion header={<h3 id="SpatialMarkers" ref={MarkerAnchor}>Spatial Markers</h3>}>
                {<div>
                  (left) These are the assignments of each spot-barcode to clusters by an
                  <br />
                </div>}
              </ToggleAccordion>
              <MarkerTable file={statusSP.data} setCompLoad={setComponentLoad} onRef={MarkerRef} />
            </Card>
            <Divider />
            <Card>
              <ToggleAccordion header={<h3 id="CellInteractions" ref={CPDBAnchor}>Cell Interactions</h3>}>
                {<div>
                  (left) These are the assignments of each spot-barcode to clusters by an
                  <br />
                  (right) Spots are colored by clustering assignment and shown in t-SNE space.
                </div>}
              </ToggleAccordion>
              <NetworkRelation spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} onRef={NetRef} />
              <RelHeat spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} onRef={CPDBRef} />
            </Card>
            {/* <div style={{color: 'black'}}>{JSON.stringify(location.state)}</div> */}
          </Content>
        </Layout>
      </Layout>
    </ConfigProvider>
  )
}
export default Browser