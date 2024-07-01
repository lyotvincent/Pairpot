import { DotChartOutlined, FileDoneOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, SelectOutlined, SlidersOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense, useRef } from 'react'
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
import { Cursor } from 'react-bootstrap-icons'
const { Content, Sider, Header } = Layout

const sideMenuItems = [
  {
    label: 'Search',
    href: '#Search',
    key: 'Search',
    icon: <FileSearchOutlined />
  },
  {
    label: 'MetaInfo',
    href: '#MetaInfo',
    key: 'MetaInfo',
    icon: <FileDoneOutlined />
  },
  {
    label: 'LassoView',
    href: '#LassoView',
    key: 'LassoView',
    icon: <SelectOutlined />
  },
  {
    label: 'Layer View',
    href: '#LayerView',
    key: 'LayerView',
    icon: <DotChartOutlined />
  },
  {
    label: 'Pair View',
    href: '#PairView',
    key: 'PairView',
    icon: <OneToOneOutlined />
  },
  {
    label: 'Spatial Markers',
    href: '#SpatialMarkers',
    key: 'SpatialMarkers',
    icon: <SlidersOutlined />,
  },
  {
    label: 'Cell Interactions',
    href: '#CellInteractions',
    key: 'CellInteractions',
    icon: <HeatMapOutlined />,
  },
]

const MetaInfo = lazy(() => (import('./utils/DatasetDetails')))
const LassoView = lazy(() => (import('./charts/ScScatter')))
const LayerView = lazy(() => (import('./charts/SpScatter')))
const PairView = lazy(() => (import('./charts/PairView')))

const Browser = () => {
  const navigate = useNavigate()
  const location = useLocation()
  const [Init, setInit] = useState(false)
  const [collapsed, setCollapsed] = useState(false)
  const [selectedKey, setSelectedKey] = useState(sideMenuItems[0].key)
  const [action, setAction] = useState(0)
  const commandRef = useRef()
  const onClickSideMenu = (e) => {

  }
  const fetchSPData = () => (
    axios.get('/api/query/sp', {
      params: {
        id: location.state?.st?.dataset_id
      },
      responseType: 'blob',
    }).then((response) => {
      let spblob = response.data
      let _spfile = new File([spblob], "sp_meta.h5ad")
      return _spfile
    })
  )

  const fetchSCData = (loc) => (
    axios.get('/api/query/sc', {
      params: {
        id: loc.state?.st?.dataset_id
      },
      responseType: 'blob',
    }).then((response) => {
      let scblob = response.data
      let _scfile = new File([scblob], "sc_sampled.h5ad")
      return _scfile
    })
  )

  const LassoRef = React.createRef()  // get the trigger function from LassoView
  const LayerRef = React.createRef()  // get the trigger function from LayerView
  const PairRef = React.createRef()   // get the trigger function from PairView
  const MarkerRef = React.createRef()   // get the trigger function from Marker Table
  const NetRef = React.createRef()   // get the trigger function from Network Table
  const CPDBRef = React.createRef()   // get the trigger function from CPDB heatmap

  const sc = fetchSCData(location)
  const sp = fetchSPData(location)

  const [scfile, setSCFile] = useState(sc)
  const [spfile, setSPFile] = useState(sp)

  const toggleAction = (command) => {
    commandRef.current = command
    setAction(action + 1)
  }

  useEffect(() => {
    console.log(location.state)
    if (Init && location.state !== null) {
      setSCFile(fetchSCData(location))
      setSPFile(fetchSPData(location))
      LassoRef.current.Trigger("Reload")  // call ToggleAnno to trigger reload in useEffect
      LayerRef.current.Trigger("Reload")
      PairRef.current.Trigger("Reload")
      MarkerRef.current.Trigger("Reload")
      NetRef.current.Trigger("Reload")
      CPDBRef.current.Trigger("Reload")
    }
    else {
      if (location.state === null) {
        axios({
          method: 'GET',
          url: '/api/example',
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
          navigate('/browse', { state: state })
        })
      } else {
        setInit(true)
      }
    }
  }, [location])

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
            <Anchor targetOffset={100} affix={false} >
              {sideMenuItems.map((item) => (
                <Menu.Item key={item.key} icon={item.icon} style={{ paddingLeft: "1.5rem" }}>
                  <Anchor.Link href={item.href} title={item.label} >
                  </Anchor.Link>
                </Menu.Item>
              ))}
            </Anchor>
          </Menu>
        </Sider>
        <Layout>
          <Content style={{ padding: '20px 20px' }}>
            <Col span={6} offset={8} id='Search'>
              <div style={{
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
            <Suspense fallback={<h1>Loading for MetaInfo...</h1>}>
              <Card id="MetaInfo">
                <MetaInfo
                  location={location}
                  navigate={navigate} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LassoView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="LassoView">Lasso View</h3>}>
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
                <LassoView spfile={spfile} scfile={scfile} location={location} onRef={LassoRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for LayerView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="LayerView">Layer View</h3>}>
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
                <LayerView id='LayerView' query={true} spfile={spfile} onRef={LayerRef} />
              </Card>
            </Suspense>
            <Divider />
            <Suspense fallback={<h1>Loading for PairView...</h1>}>
              <Card>
                <ToggleAccordion header={<h3 id="PairView">Pair View</h3>}>
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
                <PairView spfile={spfile} scfile={scfile} location={location} onRef={PairRef} />
              </Card>
            </Suspense>
            <Divider />
            <Card>
              <ToggleAccordion header={<h3 id="SpatialMarkers">Spatial Markers</h3>}>
                {<div>
                  (left) These are the assignments of each spot-barcode to clusters by an
                  <br />
                </div>}
              </ToggleAccordion>
              <MarkerTable file={spfile} onRef={MarkerRef} />
            </Card>
            <Divider />
            <Card>
              <ToggleAccordion header={<h3 id="CellInteractions">Cell Interactions</h3>}>
                {<div>
                  (left) These are the assignments of each spot-barcode to clusters by an
                  <br />
                  (right) Spots are colored by clustering assignment and shown in t-SNE space.
                </div>}
              </ToggleAccordion>
              <NetworkRelation spfile={spfile} scfile={scfile} onRef={NetRef} />
              <RelHeat spfile={spfile} scfile={scfile} onRef={CPDBRef} />
            </Card>
            {/* <div style={{color: 'black'}}>{JSON.stringify(location.state)}</div> */}
          </Content>
        </Layout>
      </Layout>
    </ConfigProvider>
  )
}
export default Browser