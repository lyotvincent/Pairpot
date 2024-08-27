import { DotChartOutlined, FileDoneOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, SelectOutlined, SlidersOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense, useRef, memo } from 'react'
import logoFig from "../assets/img/mylogo.png"
import {
  Menu,
  Button,
  Input,
  Col,
  Form,
  Steps,
  Row,
  Space,
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
import contactImg from "../assets/img/contact.png"
import { useMutation, useQuery, QueryCache } from "react-query"
const { Content, Sider, Header, Footer } = Layout
const { enterLoading } = Loading

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
    queryFn: () => axios.get('/api/example', {
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

  useEffect(() => {
    if (example.status === 'success' && typeof example.data !== 'undefined' && location.state === null) {
      setLocState(example.data)
    }
  }, [example.data])

  useEffect(() => {
    if (typeof locState !== 'undefined' &&
      location.state === null &&
      MetaRef.current !== null) {
      MetaRef.current.Trigger(example.data)
      statusSC.mutate(locState)
      statusSP.mutate(locState)
    }
  }, [locState, MetaRef.current])

  useEffect(() => {
    if (Init && location.state !== null) {
      // set MetaInfo
      MetaRef.current.Trigger(location.state)
      setLocState(location.state)
      setComponentLoad({  // reset ComponentLoad
        "LassoView": false,
        "LayerView": false,
        "PairView": false,
        "MarkerTable": false,
        "NetWorkRelation": false,
        "CPDBHeatmap": false,
      })

      // set loading status for each components
      LassoRef.current?.Tip(loadingTips[3])
      if (LassoRef.current !== null)
        enterLoading(0, LassoRef.current.Loading)

      LayerRef.current?.Tip(loadingTips[3])
      if (LayerRef.current !== null)
        enterLoading(0, LayerRef.current.Loading)

      PairRef.current?.Tip(loadingTips[3])
      if (PairRef.current !== null)
        enterLoading(1, PairRef.current.Loading)

      MarkerRef.current?.Tip(loadingTips[3])
      if (MarkerRef.current !== null) {
        enterLoading(0, MarkerRef.current.Loading)
        enterLoading(1, MarkerRef.current.Loading)
      }

      CPDBRef.current?.Tip(loadingTips[3])
      if (CPDBRef.current !== null)
        enterLoading(0, CPDBRef.current.Loading)

      NetRef.current?.Tip(loadingTips[3])
      if (NetRef.current !== null)
        enterLoading(0, NetRef.current.Loading)

      // mutation of datasets
      statusSC.mutate(location.state)
      statusSP.mutate(location.state)
    }
  }, [location, Init])

  useEffect(() => {
    if (PairRef.current !== null &&
      LassoRef.current !== null &&
      LayerRef.current !== null &&
      MetaRef.current !== null &&
      !Init
    ) {
      setInit(true)  // set Browser initializtion finished afer all components.
    }
  }, [PairRef.current, LassoRef.current, LayerRef.current, MetaRef.current])

  useEffect(() => {
    if (statusSP.status === 'success' &&
      typeof statusSP.data !== "undefined") {
      if (!componentLoad['LayerView']) {
        LayerRef.current?.Trigger("Reload")  // only reload once
      }
      if (!componentLoad['MarkerTable']) {
        MarkerRef.current?.Trigger("Reload")
      }
    }
    if (statusSP.status === 'success' &&
      statusSC.status === 'success' &&
      typeof statusSC.data !== "undefined" &&
      typeof statusSP.data !== "undefined") {
      if (!componentLoad['LassoView']) {
        LassoRef.current?.Trigger("Reload")
      }
      if (!componentLoad['PairView']) {
        PairRef.current?.Trigger("Reload")
      }
      if (!componentLoad['NetWorkRelation']) {
        NetRef.current?.Trigger("Reload")
      }
      if (!componentLoad['CPDBHeatmap']) {
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
                {/* <div>{JSON.stringify(location.state)}</div> */}
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LassoView...</h1>}>
              <Card>
                <ToggleAccordion header={
                  <h3 id="LassoView" ref={LassoAnchor}>Lasso-View
                  </h3>}>
                  {<div>
                    <b>Overview:</b> Lasso-View is a heuristic analysis that discovers extra omitted
                    cells (spots) similar to the lassoed ones in both single-cell
                    and SRT datasets based on semi-supervised machine learning.
                    <br />
                    <b>Left Chart:</b> The embeddings of single-cell or SRT data with predefined annotations
                    of spatial domains and cell types.
                    <br />
                    <b>Right Chart:</b> The embeddings of single-cell or SRT data with customized annotations according to user-selction.
                    <br />
                    <b>Usage:</b>
                    <ol>
                      <li>Click the lasso tools at TopRight, and then select the
                        cells of interests in the left chart. </li>
                      <li>Use select modes and tools to zoom, draw, and erase the
                        customized cells. </li>
                      <li>Click the `Rename` button and type in the name of your selected cells. Your annotations are in the legends (bottom of the charts).</li>
                      <li>Click the `Refine` button, choose a refiner (currently `LabelPropagation` is available), and then click `Refine` button to refine your annotations.</li>
                      <li>Click `Confirm` button to ensure your current annotation. Then you can lasso for other annotations.</li>
                      <li>Lasso your new anntations as same as step 1. You can delete your annotations through `Delete` button.</li>
                      <li>Save your annotations as a JSON file.</li>
                    </ol>
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
                <LassoView spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} meta={locState} onRef={LassoRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LayerView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="LayerView" ref={LayerAnchor}>Layer-View</h3>}>
                  {<div>
                    <b>Overview:</b> Inspired by 3D Landscape in De-spot, Pairpot develops
                    Layer-View to dynamically display gene expressions and
                    annotations of multiple slices in 3D hierarchical layouts
                    Users can explore all slices of a study in the
                    left chart, and explore a single slice precisely in the right chart.In the 3D layout, the x-axis and y-axis denote the
                    spatial coordinates of the slices, and the z-axis represents
                    different batches of these slices.
                    <br />
                    Users can rotate the axis to
                    different perspectives according to their interests. Users can
                    click ’Inverse’ to hide all annotations and strengthen specific
                    domains by clicking their legends.
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
                <LayerView id='LayerView' query={true} spfile={statusSP.data} setCompLoad={setComponentLoad} onRef={LayerRef} meta={locState} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for PairView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="PairView" ref={PairAnchor}>Pair-View</h3>}>
                  {<div>
                    <b>Overview:</b> Pair-View is another heuristic analysis that aids users to quickly
                    infer cell proportions of spots using their customized cells
                    from single-cell data in real time. Pair-
                    View API performs NNLS online based on the UCell scores
                    pre-analyzed by Pair-View, which uses the Mann-Whitney U
                    statistic to evaluate the related gene expression in single-cell
                    and SRT data.
                    <br/>
                    <b>Left Chart:</b> The embedding of single-cell data with pre-defined cell types.
                    <br/>
                    <b>Right Chart:</b> The embedding of SRT data with pre-defined spatial domains.
                    <br/>
                    <b>Usage:</b>
                    <ol>
                      <li>Click the Lasso Tools at TopRight, then select the cells of interests in the left chart.</li>
                      <li>Use Select modes and tools to zoom, draw, and erase the customized cells.</li>
                      <li>Click the `Deconv` buttons to infer cell proportions of spots using the user-selected cells. It would cost a few seconds. </li>
                      <li>Users can also refine the selected current cell types as same as the operations in LassoView.</li>
                      <li>Drag the bar at the top of charts to adjust the thresholds.</li> 
                      <li>Explore cell proportions of different batches through `Batches` in `Settings` button.</li>
                      <li>Click `Save` to download the customized cell types and their proportions.</li>
                    </ol>
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
                  <b>Top:</b> List of spatial markers in each cluster. The clusters are segmented by MENDER in SRT data, and are indentified by Leiden in single-cell data. Users can change the number of markers for each cluster.
                  <br/>
                  <b>Bottom:</b> Heatmap of spatial markers in each cluster. Rows are genes and columns are clusters. 
                  <br />
                </div>}
              </ToggleAccordion>
              <MarkerTable file={statusSP.data} setCompLoad={setComponentLoad} onRef={MarkerRef} />
            </Card>
            <Divider />
            <Card>
              <ToggleAccordion header={<h3 id="CellInteractions" ref={CPDBAnchor}>Cell Interactions</h3>}>
                {<div>
                  <b>Top:</b> Cell-Cell interaction networks. Click the switch to change the chart from single-cell to SRT data. Click the legends to rerender the chart.
                  <br />
                  <b>Bottom:</b> Ligand-Receptor pairs. Click the switch to change the chart from single-cell to SRT data. Select the cluster ID to display its related L-R pairs. 
                </div>}
              </ToggleAccordion>
              <Row>
              <Col span={12} offset={5}>
              <NetworkRelation spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} onRef={NetRef} />
              </Col>
              </Row>
              <RelHeat spfile={statusSP.data} scfile={statusSC.data} setCompLoad={setComponentLoad} onRef={CPDBRef} />
            </Card>
            {/* <div style={{color: 'black'}}>{JSON.stringify(location.state)}</div> */}
          </Content>
          <Footer>
            <Row>
              <Col span={4}>
                <Space direction='horizontal'>
                  <h3>Contact us
                    {' '}
                    <img src={contactImg} style={{ width: '20%', height: '20%' }}></img>
                  </h3>
                </Space>
                <a href="https://beian.miit.gov.cn" target='_blank'
                  style={{
                    color: 'rgba(0, 0, 0, 0.6)',
                    textDecoration: 'none',
                    fontSize: 13
                  }}>
                  津ICP备2024022781号-1
                </a>
              </Col>
              <Col offset={1}>
                <p>
                  *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                  Medicine, Nankai University, rrrzhan@nankai.edu.cn
                </p>
                <p>
                  *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@nankai.edu.cn
                </p>
              </Col>
            </Row>
          </Footer>
        </Layout>
      </Layout>
    </ConfigProvider>
  )
}
export default Browser