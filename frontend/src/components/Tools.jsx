import { useState } from "react"
import TrackVisibility from 'react-on-screen'
import '../App.css'
import { Layout, Menu, theme, Space, Button, ConfigProvider, Card, Row, Col } from 'antd'
import 'bootstrap/dist/css/bootstrap.min.css'
import ToggleAccordion from "./utils/ToggleAccordion"
//import { FileUpload } from "./utils/FileUpload"
import SpScatter from "./charts/SpScatter"
import Landscape3D2 from "./charts/Landscape3D2"
import ScScatter from "./charts/ScScatter"
import ScHeatmap from './charts/ScHeatmap'
import DatasetStatic from './utils/DatasetStatic'
import NetworkRelation from "./charts/NetworkRelation"
// import RelationHeatmap from "./charts/Relation-Heatmap"
import RelHeat from "./charts/Relation-Heatmap-test"
import SeriesGallery from "./utils/SeriesGallery"
import CellProps from "./charts/CellProps"
import CpdbCorr from "./charts/CpdbCorr"
import axios from "axios"
import { MenuUnfoldOutlined, MenuFoldOutlined, ProfileOutlined } from "@ant-design/icons"
const { Header, Content, Footer, Sider } = Layout
const sideMenuItems = [
  {
    label: 'CellPropagation',
    key: 'CellPropagation',
    icon: <ProfileOutlined />,
  },
  {
    label: 'LayerLinkage',
    key: 'LayerLinkage',
    icon: <ProfileOutlined />,
  },
  {
    label: 'LassoDeconv',
    key: 'LassoDeconv',
    icon: <ProfileOutlined />,
  },
]

export const Tools = () => {
  const [collapsed, setCollapsed] = useState(false)
  const [selectedKey, setSelectedKey] = useState('CellPropagation')
  const [visible, setVisible] = useState(true)
  const onClickSideMenu = (e) => {
    setSelectedKey(e.key)
  }
  const fetchSPData = () => (
    axios.get('/api/query/sp', {
      responseType: 'blob',
    }).then((response) => {
      let spblob = response.data
      let _spfile = new File([spblob], "sp1_meta.h5ad")
      return _spfile
    })
  )

  return (
    <ConfigProvider theme={{
      cssVar: true,
    }}>
      <Layout hasSider
        style={{
          minHeight: '100vh',
        }}>
        <Sider
          width={200}
          breakpoint="lg"
          collapsible collapsed={collapsed} onCollapse={(value) => setCollapsed(value)}>
          <Menu
            mode="inline"
            onClick={onClickSideMenu}
            defaultSelectedKeys={'sub1'}
            style={{
              width: collapsed ? 80 : 200,
              height: '100%',
              borderRight: 0,
              position: 'fixed'
            }}
            items={sideMenuItems}
          />
        </Sider>
        <Layout>

          <Content
            style={{
              padding: '1rem 1rem',
            }}>
            {visible && selectedKey === 'CellPropagation' && (
              <Card >
                <ToggleAccordion header={<h3>Cell Propagation</h3>} border={null}>
                  {<div>
                    (left) These are the assignments of each spot-barcode to clusters by an
                    automated clustering algorithm. The clustering groups together spots that
                    have similar expression profiles. In this plot, spots are colored according
                    to their cluster assignment and projected on to the tissue image. Only spots
                    under tissue are used in the clustering algorithm.
                    <br />
                    (right) Spots are colored by clustering assignment and shown in t-SNE space.
                    The axes correspond to the 2-dimensional embedding produced by the t-SNE
                    algorithm. In this space, pairs of spots that are close to each other have
                    more similar gene expression profiles than spots that are distant from each other.
                  </div>}
                </ToggleAccordion>
                <Row>
                  <Col><ScScatter /></Col>
                  <Col ><ScHeatmap /></Col>
                  <Col ><CpdbCorr /></Col>
                  {/* <Col><NetworkRelation/></Col> */}
                </Row>

                <Row>
                  <Col><NetworkRelation file={fetchSPData()} /></Col>
                </Row>

                <Row>
                  <Col><RelHeat file={fetchSPData()} /></Col>
                </Row>
              


                <br />
                <br />
              </Card>
            )}

            <br />
            <Card >
              <ToggleAccordion header={<h3>LayerLinkages</h3>} border={null}>
                {<div>
                  (left) Total UMI counts for each spot overlayed on the tissue
                  image. Spots with greater UMI counts likely have higher RNA
                  content than spots with fewer UMI counts.
                  <br />
                  (right) Total UMI counts for spots displayed by a 2-dimensional
                  embedding produced by the t-SNE algorithm. In this space, pairs of
                  spots that are close to each other have more similar gene
                  expression profiles than spots that are distant from each other.
                </div>}
              </ToggleAccordion>
              <Row>
                <Col><SpScatter /></Col>
              </Row>
            </Card>

            <br />
            <Card >
              <ToggleAccordion header={<h3>Cell Propagation</h3>} border={null}>
                {<div>
                  (left) These are the assignments of each spot-barcode to clusters by an
                  automated clustering algorithm. The clustering groups together spots that
                  have similar expression profiles. In this plot, spots are colored according
                  to their cluster assignment and projected on to the tissue image. Only spots
                  under tissue are used in the clustering algorithm.
                  <br />
                  (right) Spots are colored by clustering assignment and shown in t-SNE space.
                  The axes correspond to the 2-dimensional embedding produced by the t-SNE
                  algorithm. In this space, pairs of spots that are close to each other have
                  more similar gene expression profiles than spots that are distant from each other.
                </div>}
              </ToggleAccordion>
              <Row>
                <Col><CellProps /></Col>
              </Row>
            </Card>

            <br />

            <Card >
              <ToggleAccordion header={<h3>3D Landscapes</h3>} border={null}>
                {<div>
                  The generated landscapes curve the outlines of SCSPs in tissues, and display the
                  overlapped SCSPs in layered structures. Despot uses the Alpha-shape algorithm to curve
                  the outlines in tissues and uses dots with colors to represent the layers of SCSPs.
                  For each SCSP layer, Despot randomly select dots of 10 percent and create waterfall
                  lines to link the dots and outlined regions together.
                </div>}
              </ToggleAccordion>
              <Row>
                <Landscape3D2 />
              </Row>
            </Card>
            <br />
            <Card >
              <ToggleAccordion header={<h3>Data Selection & Submittion</h3>} border={null}>
                {<div>
                  The Selection module supports .h5spt, .h5ad, .mtx or .zip file with 10X Visium format.
                  Hence, .mtx needs to have genes in rows and barcodes(cells) in columns.
                </div>}
              </ToggleAccordion>
              <h4 align="center"> Series Gallery</h4>
              <SeriesGallery />
            </Card>
          </Content>
        </Layout>
      </Layout>
    </ConfigProvider>
  )
}