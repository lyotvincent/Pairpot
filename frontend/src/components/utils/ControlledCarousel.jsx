import React, { useState } from 'react'
import Carousel from 'react-bootstrap/Carousel'
import DotPlot from '../charts/DotPlot'
import BoxPlot from '../charts/BoxPlot'
import ScatterPlot from '../charts/ScatterPlot'

export const ControlledCarousel = () => {
  const [index, setIndex] = useState(0)

  const handleSelect = (selectedIndex, e) => {
    setIndex(selectedIndex)
  }

  return (
    <Carousel activeIndex={index} onSelect={handleSelect}>
      <Carousel.Item>
        <DotPlot />
        <Carousel.Caption>
          <h3>Spatial Clustering</h3>
          <p>Nulla vitae elit libero, a pharetra augue mollis interdum.</p>
        </Carousel.Caption>
      </Carousel.Item>
      <Carousel.Item>
        <BoxPlot />
        <Carousel.Caption>
          <h3>Spatial Variable Genes</h3>
          <p>Lorem ipsum dolor sit amet, consectetur adipiscing elit.</p>
        </Carousel.Caption>
      </Carousel.Item>
      <Carousel.Item>
        <ScatterPlot />
        <Carousel.Caption>
          <h3>Deconvolution</h3>
          <p>
            Praesent commodo cursus magna, vel scelerisque nisl consectetur.
          </p>
        </Carousel.Caption>
      </Carousel.Item>
    </Carousel>
  )
}
