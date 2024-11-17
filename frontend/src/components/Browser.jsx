import { DotChartOutlined, FileDoneOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, SelectOutlined, SlidersOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense, useRef, memo } from 'react'
import { ZipReader, BlobReader, BlobWriter } from "@zip.js/zip.js"
import logoFig from "../assets/img/mylogo.png"
import githubLogo from "../assets/img/github.svg"
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
  Segmented,
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
  const testFlag = true
  const navigate = useNavigate()
  const location = useLocation()
  const [Init, setInit] = useState(false)
  const [collapsed, setCollapsed] = useState(false)
  const [locState, setLocState] = useState(false)
  const [action, setAction] = useState(0)
  const commandRef = useRef()
  const stProgress = useRef()
  const scProgress = useRef()
  const onClickSideMenu = (e) => {

  }
  const [myVal, setMyval] = useState(location.state != null ? location.state.myVal : 'null')
  // const myVal = () => {
  //   if(location.state == null)
  //     return 'null'
  //   if(location.state.myVal == 'LassoView' ||
  //      location.state.myVal == "Layer View" ||
  //      location.state.myVal == "Pair View" || 
  //      location.state.myVal == "Cell Interactions"
  //   ){
  //     const val = location.state.myVal
  //     // 使用 navigate 只更新 state，不改变路径
  //     // navigate('.', { replace: true, state: null });
  //     return val
  //   }
  // }

  const [modelKey, setModelKey] = useState(true)

  const fetchSPData = (state) => (
    axios.get('/api/query/spgz', {
      params: {
        id: state.st?.dataset_id
      },
      responseType: 'blob',
      onDownloadProgress: (progressEvent) => {
        stProgress.current = Math.round(progressEvent.loaded / progressEvent.total * 100)
      }
    }).then((response) => {
      let spblob = response.data
      let zipFileReader = new BlobReader(spblob)
      let zipReader = new ZipReader(zipFileReader)
      let firstEntry = zipReader.getEntries().then((entry) => { return (entry.shift()) })
      return (
        firstEntry.then((entry) => entry.getData(new BlobWriter()))
          .then((_blob) => {
            return (new File([_blob], "sp_meta.h5ad"))
          })
      )
    }).then((_scfile) => {
      return (_scfile)
    })
  )

  const fetchSCData = (state) => (
    axios.get('/api/query/scgz', {
      params: {
        id: state.st?.dataset_id
      },
      responseType: 'blob',
      onDownloadProgress: (progressEvent) => {
        scProgress.current = Math.round(progressEvent.loaded / progressEvent.total * 100)
      }
    }).then((response) => {
      let scblob = response.data
      let zipFileReader = new BlobReader(scblob)
      let zipReader = new ZipReader(zipFileReader)
      let firstEntry = zipReader.getEntries().then((entry) => { return (entry.shift()) })
      return (
        firstEntry.then((entry) => entry.getData(new BlobWriter()))
          .then((_blob) => {
            return (new File([_blob], "sc_meta.h5ad"))
          })
      )
    }).then((_scfile) => {
      return (_scfile)
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
      // console.log("get response")
      // console.log(response)
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
  const heatAnchor = React.createRef() // location of heat-net
  const sideMenuItems = [
    {
      label: 'MetaInfo',
      ref: MetaAnchor,
      key: 'MetaInfo',
      icon: <FileDoneOutlined />
    },
    {
      label: 'Lasso-View',
      ref: LassoAnchor,
      key: 'LassoView',
      icon: <SelectOutlined />
    },
    {
      label: 'Layer-View',
      ref: LayerAnchor,
      key: 'LayerView',
      icon: <DotChartOutlined />
    },
    {
      label: 'Pair-View',
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
      label: 'Interactions',
      ref: CPDBAnchor,
      key: 'Interactions',
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
    // console.log(location.state)
    // console.log(example)
    if (example.status === 'success' && typeof example.data !== 'undefined' &&
      (location.state === null || myVal !== 'null')) {
      // console.log("setting example")
      setLocState(example.data)
      // console.log("example finished")
    }
  }, [example.data])

  useEffect(() => {
    if (typeof locState !== 'undefined' && locState !== null &&
      (location.state === null || myVal !== 'null') &&
      MetaRef.current !== null) {
      MetaRef.current.Trigger(locState)
      statusSC.mutate(locState)
      statusSP.mutate(locState)
    }
  }, [locState, MetaRef.current])

  useEffect(() => {
    if (Init && location.state !== null) {
      setLocState(location.state)
    }
  }, [location])

  useEffect(() => {
    // console.log("Init:", Init)
    if (Init && locState !== null && myVal === 'null') {
      // set MetaInfo
      // console.log("Setting MetaInfo")
      MetaRef.current.Trigger(locState)
      setLocState(locState)
      setComponentLoad({  // reset ComponentLoad
        "LassoView": false,
        "LayerView": false,
        "PairView": false,
        "MarkerTable": false,
        "NetWorkRelation": false,
        "CPDBHeatmap": false,
      })

      // set loading status for each components
      LassoRef.current?.Tip(<div>{loadingTips[3]}</div>)
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

      // console.log("netref")
      console.log(NetRef)
      NetRef.current?.Tip(loadingTips[3])
      if (NetRef.current !== null)
        enterLoading(0, NetRef.current.Loading)

      // mutation of datasets
      // console.log(locState)
      statusSC.mutate(locState)
      statusSP.mutate(locState)
    }
  }, [locState, Init])

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
    if (!modelKey) {
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
    }
    else {
      if (statusSP.status === 'success' &&
        statusSC.status === 'success' &&
        typeof statusSC.data !== "undefined" &&
        typeof statusSP.data !== "undefined") {
        if (selectedKey === sideMenuItems[1].key &&
          !componentLoad['LassoView']
        ) {  // LassoView
          LassoRef.current?.Trigger("Reload")
        }
        if (selectedKey === sideMenuItems[2].key &&
          !componentLoad['LayerView']
        ) {  // LayerView
          LayerRef.current?.Trigger("Reload")
        }
        if (selectedKey === sideMenuItems[3].key &&
          !componentLoad['PairView']
        ) {  // PairView
          PairRef.current?.Trigger("Reload")
        }
        if (selectedKey === sideMenuItems[4].key &&
          !componentLoad['MarkerTable']
        ) {  // MarkerTable
          MarkerRef.current?.Trigger("Reload")
        }
        if (selectedKey === sideMenuItems[5].key &&
          !componentLoad['NetWorkRelation']
        ) {  // NetWorkRelation
          NetRef.current?.Trigger("Reload")
        }
        if (selectedKey === sideMenuItems[5].key &&
          !componentLoad['CPDBHeatmap']
        ) {  // CPDBHeatmap
          CPDBRef.current?.Trigger("Reload")
        }
      }
    }
  }, [statusSP.variables, statusSC.variables, PairRef])

  useEffect(() => {
    // console.log(myVal)
    if (
      LassoAnchor.current != null &&
      LayerAnchor.current != null &&
      PairAnchor.current != null &&
      CPDBAnchor.current != null &&
      SearchAnchor.current != null
      || true // only need true
    ) {
      if (myVal != 'null' && myVal != 'undefined') {
        switch (myVal) {
          case 'Lasso-View':
            // setMyValItem(sideMenuItems[2]);
            setSelectedKey(sideMenuItems[1].key)
            // sideMenuItems[2].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
            break
          case 'Layer-View':
            // setMyValItem(sideMenuItems[3]);
            // console.log(myValItem)
            setSelectedKey(sideMenuItems[2].key)
            // sideMenuItems[3].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
            break
          case 'Pair-View':
            // setMyValItem(sideMenuItems[4]);
            setSelectedKey(sideMenuItems[3].key)
            // sideMenuItems[4].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
            break
          case 'Interactions':
            // setMyValItem(sideMenuItems[6]);
            setSelectedKey(sideMenuItems[5].key)
            // sideMenuItems[6].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
            break
          case 'Spatial Markers':
            setSelectedKey(sideMenuItems[4].key)
            break;
          default:
            // setMyValItem(sideMenuItems[0]);
            setSelectedKey(sideMenuItems[0].key)
            // sideMenuItems[0].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
            break
        }

      }
    }

    // // if --> async
    // // conditions
    // const checkCondition = () => {
    //   return (
    //     PairRef.current !== null &&
    //     LassoRef.current !== null &&
    //     LayerRef.current !== null &&
    //     MetaRef.current !== null
    //   );
    // };

    // const waitForCondition = async () => {
    //   while (!checkCondition()) {
    //     console.log("Waiting for refs to be ready...");
    //     await new Promise((resolve) => setTimeout(resolve, 100)); // 每隔 100ms 检查一次
    //   }
    //   // 条件满足后执行后续操作
    //   console.log("All refs are ready!");
    //   // 这里写你的逻辑，当 if 条件满足时执行

    //   if( myVal != 'null' && myVal != 'undefined'){
    //     console.log("here")
    //     switch(myVal){
    //       case 'LassoView':
    //         // setMyValItem(sideMenuItems[2]);
    //         setSelectedKey(sideMenuItems[2].key)
    //         sideMenuItems[2].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
    //         break;
    //       case 'Layer View':
    //         // setMyValItem(sideMenuItems[3]);
    //         // console.log(myValItem)
    //         setSelectedKey(sideMenuItems[3].key)
    //         sideMenuItems[3].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
    //         break;
    //       case 'Pair View':
    //         // setMyValItem(sideMenuItems[4]);
    //         setSelectedKey(sideMenuItems[4].key)
    //         sideMenuItems[4].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
    //         break;
    //       case 'Cell Interactions' :
    //         // setMyValItem(sideMenuItems[6]);
    //         setSelectedKey(sideMenuItems[6].key)
    //         sideMenuItems[6].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
    //         break;
    //       default:
    //         // setMyValItem(sideMenuItems[0]);
    //         setSelectedKey(sideMenuItems[0].key)
    //         sideMenuItems[0].ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
    //         break;
    //     }

    //   }

    // };

    // waitForCondition(); // 启动异步轮询


  }, [myVal, Init])

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
            {/* <Switch
            checkedChildren="card"
            unCheckedChildren="roll" */}
            <Segmented block
              options={[
                { label: 'card', value: true },
                { label: 'roll', value: false },
              ]}
              value={modelKey}
              // style={{
              //   zIndex: 1000,
              //   left: '50px',
              //   top: '680px',
              //   position: 'fixed'
              // }}
              onChange={(value) => {
                setModelKey(value)
              }}
            />
            {sideMenuItems.map((item) => (
              <Menu.Item
                key={item.key}
                icon={item.icon}
                onClick={(e) => {
                  setSelectedKey(e.key)
                  item.ref.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                  console.log(e.key)
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
            <Col span={18} offset={3}>
              <Search setLocState={setLocState}
                setCompLoad={setComponentLoad}
                setSelectedKey={setSelectedKey} />
              {/* <div>{JSON.stringify(myVal)}</div> */}
            </Col>
            <Divider />
            {/* <div>compLoad:{JSON.stringify(componentLoad)}</div> */}
            <Suspense fallback={<h1>Loading for MetaInfo...</h1>}>
              <Card style={{ display: selectedKey == "MetaInfo" || !modelKey ? 'block' : 'none' }}
                id="MetaInfo" ref={MetaAnchor}>
                <MetaInfo onRef={MetaRef} />
                {/* <div>{JSON.stringify(location.state)}</div> */}
              </Card>
            </Suspense>
            <Divider style={{ display: !modelKey ? 'block' : 'none' }} />
            <Suspense fallback={<h1>Loading for LassoView...</h1>}>
              <Card style={{ display: selectedKey == "LassoView" || !modelKey ? 'block' : 'none' }}>
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
                <LassoView
                  spfile={statusSP.data}
                  scfile={statusSC.data}
                  setCompLoad={setComponentLoad}
                  meta={locState}
                  onRef={LassoRef}
                  progress={(stProgress.current + scProgress.current) / 2} />
              </Card>
              {/* )} */}
            </Suspense>
            <Divider style={{ display: !modelKey ? 'block' : 'none' }} />
            <Suspense fallback={<h1>Loading for LayerView...</h1>}>
              {/* { selectedKey == "LayerView" &&( */}
              <Card style={{ display: selectedKey == "LayerView" || !modelKey ? 'block' : 'none' }}>
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
                <LayerView
                  id='LayerView'
                  query={true}
                  spfile={statusSP.data}
                  setCompLoad={setComponentLoad}
                  onRef={LayerRef}
                  meta={locState}
                  progress={(stProgress.current + scProgress.current) / 2}
                />
              </Card>
              {/* )} */}
            </Suspense>
            <Divider style={{ display: !modelKey ? 'block' : 'none' }} />
            <Suspense fallback={<h1>Loading for PairView...</h1>}>
              {/* {selectedKey=="PairView" &&( */}
              <Card style={{ display: selectedKey == "PairView" || !modelKey ? 'block' : 'none' }}>
                <ToggleAccordion header={<h3 id="PairView" ref={PairAnchor}>Pair-View</h3>}>
                  {<div>
                    <b>Overview:</b> Pair-View is another heuristic analysis that aids users to quickly
                    infer cell proportions of spots using their customized cells
                    from single-cell data in real time. Pair-
                    View API performs NNLS online based on the UCell scores
                    pre-analyzed by Pair-View, which uses the Mann-Whitney U
                    statistic to evaluate the related gene expression in single-cell
                    and SRT data.
                    <br />
                    <b>Left Chart:</b> The embedding of single-cell data with pre-defined cell types.
                    <br />
                    <b>Right Chart:</b> The embedding of SRT data with pre-defined spatial domains.
                    <br />
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
                <PairView
                  spfile={statusSP.data}
                  scfile={statusSC.data}
                  setCompLoad={setComponentLoad}
                  meta={locState}
                  onRef={PairRef}
                  progress={(stProgress.current + scProgress.current) / 2}
                />
              </Card>
            </Suspense>
            <Divider style={{ display: !modelKey ? 'block' : 'none' }} />
            {/* { (selectedKey == "SpatialMarkers")&& ( */}
            <Card style={{ display: selectedKey == "SpatialMarkers" || !modelKey ? 'block' : 'none' }}>
              <ToggleAccordion header={<h3 id="SpatialMarkers" ref={MarkerAnchor}>Spatial Markers</h3>}>
                {<div>
                  <b>Top:</b> List of spatial markers in each cluster. The clusters are segmented by MENDER in SRT data, and are indentified by Leiden in single-cell data. Users can change the number of markers for each cluster.
                  <br />
                  <b>Bottom:</b> Heatmap of spatial markers in each cluster. Rows are genes and columns are clusters.
                  <br />
                </div>}
              </ToggleAccordion>
              <MarkerTable
                file={statusSP.data}
                setCompLoad={setComponentLoad}
                onRef={MarkerRef}
                progress={(stProgress.current + scProgress.current) / 2}
              />
            </Card>
            {/* )} */}
            <Divider style={{ display: !modelKey ? 'block' : 'none' }} />
            {/* { (selectedKey == "CellInteractions") && ( */}
            <Card style={{ display: selectedKey == "Interactions" || !modelKey ? 'block' : 'none' }}>
              <ToggleAccordion header={<h3 id="Interactions" ref={CPDBAnchor}>Cell/Domain Interactions</h3>}>
                {<div>
                  <b>Top:</b> Cell-Cell interaction networks. Click the switch to change the chart from single-cell to SRT data. Click the legends to rerender the chart.
                  <br />
                  <b>Bottom:</b> Ligand-Receptor pairs. Click the switch to change the chart from single-cell to SRT data. Select the cluster ID to display its related L-R pairs.
                </div>}
              </ToggleAccordion>
              <Row>
              <Button
                  type="primary"
                  size='middle'
                  onClick={() => {
                    CPDBAnchor.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                    setTimeout(() => {
                      NetRef.current.Tour(true)
                    }, 600)
                  }}>
                  {"Begin Tour >>>"}
                </Button>
              </Row>
              <Row>
                <Col span={12} offset={5}>
                  <NetworkRelation
                    spfile={statusSP.data}
                    scfile={statusSC.data}
                    setCompLoad={setComponentLoad}
                    onRef={NetRef}
                    progress={(stProgress.current + scProgress.current) / 2}
                  />
                </Col>
              </Row>
              <Row>
              <Button
                  ref={heatAnchor}
                  type="primary"
                  size='middle'
                  onClick={() => {
                    heatAnchor.current.scrollIntoView({ behavior: "smooth", block: 'start' })
                    setTimeout(() => {
                      CPDBRef.current.Tour(true)
                    }, 600)
                  }}>
                  {"Begin Tour >>>"}
                </Button>
              </Row>
              <RelHeat
                spfile={statusSP.data}
                scfile={statusSC.data}
                setCompLoad={setComponentLoad}
                onRef={CPDBRef}
                progress={(stProgress.current + scProgress.current) / 2}
              />
            </Card>
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
              <Col>
              <a href="https://github.com/lyotvincent/Pairpot">
                  <img src={githubLogo} alt="GitHub Logo" />
              </a>
              </Col>
              <Col offset={1}>
                <p>
                  *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                  Medicine, Nankai University, rrrzhan@mail.nankai.edu.cn
                </p>
                <p>
                  *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@mail.nankai.edu.cn
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