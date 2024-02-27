import { useState } from "react"
import { Container, Row, Col, Card } from "react-bootstrap"
import TrackVisibility from 'react-on-screen'
import { Button, Stack, Badge, Form } from 'react-bootstrap'
import 'bootstrap/dist/css/bootstrap.min.css'
import HeatMap from './charts/HeatMap'
import BoxPlot from './charts/BoxPlot'
import DotPlot from "./charts/DotPlot"
import { LoadingButton } from "./utils/LoadingButton"
import { SwitchButton } from "./utils/SwitchButton"
import { ControlledCarousel } from "./utils/ControlledCarousel"
import ToggleAccordion from "./utils/ToggleAccordion"
import HistPlot from "./charts/HistPlot"
import MyCharts from "./charts/Mycharts"
import { FileUpload } from "./utils/FileUpload"
import Landscape3D from "./charts/Landscape3D"
import SpScatter from "././charts/SpScatter"
import Image3D from "./charts/Image3D"
import Landscape3D2 from "./charts/Landscape3D2"
import SpatialPlot from "./charts/SpatialPlot"
import ScScatter from "./charts/ScScatter"
import ScHeatmap from './charts/ScHeatmap'
import Drag from "./charts/Drag"
import CellProps from "./charts/CellProps"

export const Browser = () => {
  return (
    <>
      <div>
        <h2>
          Spatial Dynamic Annotation and Visualization{' '}
          <Badge bg="secondary" as="Button">
            New
          </Badge>
          <p>this is A Browser for Spatial Dynamic Annotation and Visualization</p>
        </h2>
      </div>
      <Container fluid='true'>
        <Card bg="secret" border="secondary">
          <ToggleAccordion header={<h3>Spatial Clustering</h3>} border={null}>
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
            <Col><SpScatter /></Col>
          </Row>
        </Card>

        <br />
        <Card bg="secret" border="secondary">
          <ToggleAccordion header={<h3>UMI Detection & Decontamination</h3>} border={null}>
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
            <Col ><ScHeatmap /></Col>
          </Row>
        </Card>

        <br />
        <Card bg="secret" border="secondary">
          <ToggleAccordion header={<h3>Spatial Clustering</h3>} border={null}>
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

        <Card bg="secret" border="secondary">
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
        <Card bg="secret" border="secondary">
          <ToggleAccordion header={<h3>Data Selection & Submittion</h3>} border={null}>
            {<div>
              The Selection module supports .h5spt, .h5ad, .mtx or .zip file with 10X Visium format.
              Hence, .mtx needs to have genes in rows and barcodes(cells) in columns.
            </div>}
          </ToggleAccordion>
          <FileUpload />
        </Card>
      </Container>
    </>
  )
}