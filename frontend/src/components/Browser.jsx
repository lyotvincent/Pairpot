import { DotChartOutlined, FileSearchOutlined, FileTextOutlined, HeatMapOutlined, LinkOutlined, OneToOneOutlined, UnorderedListOutlined } from '@ant-design/icons'
import React, { useState, useEffect, lazy, Suspense } from 'react'
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
} from 'antd'
import ToggleAccordion from './utils/ToggleAccordion'
import { Link, useLocation, useNavigate } from 'react-router-dom'
import MarkerTable from './charts/MarkerTable'
import axios from 'axios'
import NetworkRelation from './charts/NetworkRelation'
import RelHeat from './charts/Relation-Heatmap-test'
const { Content, Sider } = Layout

const sideMenuItems = [
  {
    label: 'MetaInfo',
    href: '#MetaInfo',
    key: 'MetaInfo',
    icon: <FileSearchOutlined />
  },
  {
    label: 'LassoView',
    href: '#LassoView',
    key: 'LassoView',
    icon: <DotChartOutlined />
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
    icon: <HeatMapOutlined />,
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
  const [collapsed, setCollapsed] = useState(false)
  const [selectedKey, setSelectedKey] = useState(sideMenuItems[0].key)
  const onClickSideMenu = (e) => {

  }
  const fetchSPData = () => (
    axios.get('http://localhost:5522/query/sp', {
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

  const fetchSCData = () => (
    axios.get('http://localhost:5522/query/sc', {
      params: {
        id: location.state?.st?.dataset_id
      },
      responseType: 'blob',
    }).then((response) => {
      let scblob = response.data
      let _scfile = new File([scblob], "sc_sampled.h5ad")
      return _scfile
    })
  )

  const sc = fetchSCData()
  const sp = fetchSPData()
  useEffect(() => {
  }, [])

  return (
    <ConfigProvider
      theme={{
        components: {
          Anchor: {
            linkPaddingInlineStart: 0
          },
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
              width: collapsed ? 80 : 200,
              height: '100%',
              borderRight: 0,
              position: 'fixed'
            }}
          >
            <Anchor targetOffset={90}>
              {sideMenuItems.map((item) => (
                <Menu.Item key={item.key} icon={item.icon}>
                  <Anchor.Link href={item.href} title={item.label}>
                  </Anchor.Link>
                </Menu.Item>
              ))}
            </Anchor>
          </Menu>
        </Sider>
        <Content style={{ padding: '20px 20px' }}>
          <Suspense fallback={<h1>Loading...</h1>}>
            <div id="MetaInfo" >
              <MetaInfo
                location={location}
                navigate={navigate} />
            </div>
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
              <LassoView spfile={sp} scfile={sc} location={location} />
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
              <LayerView id='LayerView' query={true} spfile={sp} />
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
              <PairView spfile={sp} scfile={sc} location={location} />
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
            <MarkerTable file={sp} />
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
            <NetworkRelation spfile={sp} scfile={sc} />
            <RelHeat spfile={sp} scfile={sc} />
          </Card>
          {/* <div style={{color: 'black'}}>{JSON.stringify(location.state)}</div> */}
        </Content>
      </Layout>
    </ConfigProvider>
  )
}
export default Browser