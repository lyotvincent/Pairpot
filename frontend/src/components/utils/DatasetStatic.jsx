import { SearchOutlined } from '@ant-design/icons'
import * as echarts from 'echarts'
import React, { useEffect, useRef, useState } from 'react'
import Highlighter from 'react-highlight-words'
import { Button, Input, Space, Table, Pagination } from 'antd'
import PropTypes from 'prop-types'
import { TooltipComponent, ToolboxComponent } from 'echarts/components'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
} from 'echarts/components'
import { PieChart } from 'echarts/charts'
import { UniversalTransition, LabelLayout } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import axios from 'axios'
echarts.use([
  ToolboxComponent,
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  PieChart,
  LegendComponent,
  LabelLayout,
  CanvasRenderer,
  UniversalTransition,
])

const DatasetStatic = ({ src, col, height, width, margin }) => {
  const [techSrc, setTechSrc] = useState([])
  const [techCol, setTechCol] = useState([])
  const chartRef = useRef(null) // get current DOM container

  useEffect(() => {
    if (techSrc.length === 0) {
      axios({
        method: 'GET',
        url: 'http://localhost:5522/strategies',
      }).then((response) => {
        setTechSrc(response.data.data)
        setTechCol(response.data.attributes)
        console.log(techCol)
      })
    }
    var myChart = echarts.init(chartRef.current) //init the echart container

    var straDict = {} // record the strategy of each technology
    var straIndex = techCol.indexOf('strategy')
    for (var i = 0; i < techSrc.length; i++) {
      var stra = techSrc[i][1]
      straDict[stra] = techSrc[i][straIndex]
    }

    var attrStack = {}
    for (var i = 0; i < col.length; i++) {
      var colname = col[i]
      attrStack[colname] = i
    }

    // record which technology is used in each dataset
    var techData = src.map((item) => {
      return item[attrStack.technologies]
    })

    const noSpatialTech = ['scRNA', 'snRNA', 'RNAseq', 'scATAC']
    var techProp = {} // proportions of technologies
    var straProp = {} // proportions of strategies used by technologies
    for (var item of techData) {
      if (item !== null) {
        var _items = item.split(';')
        for (var i of _items) {
          if (!noSpatialTech.includes(i)) {
            techProp[i] = techProp[i] ? techProp[i] + 1 : 1
            straProp[straDict[i]] = straProp[straDict[i]]
              ? straProp[straDict[i]] + 1
              : 1
          }
        }
      }
    }

    var techSource = []
    for (var key of Object.keys(techProp)) {
      techSource.push({ value: techProp[key], name: key })
    }
    techSource.sort((a, b) => {
      return b.value - a.value
    })

    var straSource = []
    for (var key of Object.keys(straProp)) {
      straSource.push({ value: straProp[key], name: key })
    }

    // record the species of dataset
    var speData = src.map((item) => {
      return item[attrStack.species]
    })

    var speProp = {} // proportions of species
    const mainSpecies = ['Homo sapiens', 'Mus musculus']
    for (var item of speData) {
      if (item !== null) {
        var _items = item.split(';')
        for (var i of _items) {
          var _i = mainSpecies.includes(i) ? i : 'Others'
          speProp[_i] = speProp[_i] ? speProp[_i] + 1 : 1
        }
      }
    }
    var speSource = []
    for (var key of Object.keys(speProp)) {
      speSource.push({ value: speProp[key], name: key })
    }

    // record the tissues of dataset
    var tisData = src.map((item) => {
      return item[attrStack.tissues]
    })
    var tisProp = {} // proportions of tissues
    for (var item of tisData) {
      if (item !== null) {
        var _items = item.split(';')
        for (var i of _items) {
          tisProp[i] = tisProp[i] ? tisProp[i] + 1 : 1
        }
      }
    }

    var tisSource = []
    for (var key of Object.keys(tisProp)) {
      tisSource.push({ value: tisProp[key], name: key })
    }
    tisSource.sort((a, b) => {
      return b.value - a.value
    })

    // first plot the proportion of technologies
    myChart.setOption({
      title: {
        text: 'Proportion of Datasets',
        top: 'top',
        left: 'center',
      },
      tooltip: {
        trigger: 'item',
        formatter: '<b>{a}</b> <br/> {b}:  {c}({d}%)',
      },
      series: [
        {
          name: 'Strategy',
          type: 'pie',
          center: ['28%', '35%'],
          radius: '20%',
          label: {
            show: true,
            position: 'inner',
            formatter: '{b}\n {d}%',
            fontSize: 14,
          },
          emphasis: {
            label: {
              show: true,
              fontSize: 16,
              fontWeight: 'bold',
            },
          },
          data: straSource,
        },
        {
          name: 'Technology',
          type: 'pie',
          center: ['28%', '35%'],
          radius: ['25%', '40%'],
          roseType: false,
          minAngle: 3,
          minShowLabelAngle: 4,
          itemStyle: {
            borderRadius: 3,
            borderWidth: 1,
            borderJoin: 'round',
          },
          label: {
            show: true,
            fontSize: 14,
            position: 'outside',
          },
          labelLine: {
            length: 35,
            length2: 5,
          },
          emphasis: {
            label: {
              show: true,
              fontSize: 16,
              fontWeight: 'bold',
            },
          },
          data: techSource,
        },
        {
          name: 'Species',
          type: 'pie',
          center: ['70%', '35%'],
          radius: '20%',
          itemStyle: {
            borderRadius: 3,
            borderWidth: 1,
            borderJoin: 'round',
          },
          label: {
            show: true,
            position: 'inner',
            formatter: '{b}: {d}%',
            fontSize: 14,
          },
          emphasis: {
            label: {
              show: true,
              fontSize: 16,
              fontWeight: 'bold',
            },
          },
          data: speSource,
        },
        {
          name: 'Tissues',
          type: 'pie',
          center: ['70%', '35%'],
          radius: ['25%', '40%'],
          roseType: false,
          minAngle: 3,
          minShowLabelAngle: 4,
          itemStyle: {
            borderRadius: 3,
            borderWidth: 1,
            borderJoin: 'round',
          },
          label: {
            show: true,
            position: 'outside',
            fontSize: 14,
          },
          labelLine: {
            show: true,
            length: 50,
            length2: 5,
            smooth: true,
          },
          emphasis: {
            label: {
              show: true,
              fontSize: 16,
              fontWeight: 'bold',
            },
          },
          data: tisSource,
        },
      ],
    })
  }, [src, col, techSrc, techCol])

  return (
    <div>
      <div
        ref={chartRef}
        className="chart"
        //the target DOM container needs height and width
        style={{ height: height, width: width, margin: margin }}></div>
    </div>
  )
}

DatasetStatic.defaultProps = {
  height: '35rem',
  width: '70rem',
  margin: '0rem',
}

DatasetStatic.propTypes = {
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default DatasetStatic
