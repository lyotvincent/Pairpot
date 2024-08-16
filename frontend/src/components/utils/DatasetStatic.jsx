import { SearchOutlined } from '@ant-design/icons'
import * as echarts from 'echarts'
import React, { useEffect, useRef, useState } from 'react'
import Highlighter from 'react-highlight-words'
import { Row, Col } from 'antd'
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
import { useQuery } from 'react-query'
import Js2WordCloud from 'js2wordcloud'
import ToggleAccordion from './ToggleAccordion';
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
var flag = 1
const DatasetStatic = ({ src, col, height, width, margin }) => {
  const [techSrc, setTechSrc] = useState([])
  const [techCol, setTechCol] = useState([])
  const [Init, setInit] = useState(false)
  // const [dataSrc, setDataSrc] = useState([])
  const chartRef = useRef(null) // get current DOM container
  const wordCloudRef = useRef(null)

  const wordCloudResponse = useQuery({  // change axios to reactQuery
    queryKey: ['wordcloud'],
    queryFn: () => axios.get('/api/global_words').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })


  useEffect(() => {
    if (wordCloudResponse.status === 'success' && typeof wordCloudResponse.data !== 'undefined') {
      var dataList = Object.entries(wordCloudResponse.data)
      var wc = new Js2WordCloud(wordCloudRef.current);
      wc.setOption({
        tooltip: {
          show: true,
          backgroundColor: 'rgba(0, 0, 0, 0.701961)',
          formatter: function (item) {
            if (item[1] > 12) {
              document.querySelector('.__wc_tooltip__').style.backgroundColor = 'rgba(0, 0, 0, 0.701961)';
              return item[0] + ':' + item[2];
            } else {
              document.querySelector('.__wc_tooltip__').style.backgroundColor = 'transparent';
              return '';
            }
          }
        },
        list: dataList,
        color: function () {
          const colors = ['#FF5733', '#4CAF50', '#5733FF', '#F1C40F', '#E67E22'];
          return colors[Math.floor(Math.random() * colors.length)];
        },
        fontSizeFactor: 1,                                    // 当词云值相差太大，可设置此值进字体行大小微调，默认0.1
        maxFontSize: 80,                                        // 最大fontSize，用来控制weightFactor，默认60
        minFontSize: 7,                                        // 最小fontSize，用来控制weightFactor，默认12
        tooltip: {
          show: true,                                         // 默认：false
          backgroundColor: 'rgba(0, 0, 0, 0.701961)',         // 默认：'rgba(0, 0, 0, 0.701961)'
          formatter: function (item) {                         // 数据格式化函数，item为list的一项
          }
        },
        noDataLoadingOption: {                                  // 无数据提示。
          backgroundColor: '#eee',
          text: 'data loading',
          textStyle: {
            color: '#888',
            fontSize: 14
          }
        }
      })
    }
  }, [wordCloudResponse])

  const strategyResponse = useQuery({
    queryKey: ['strategy'],
    queryFn: () => axios.get('/api/strategies').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })
  useEffect(() => {
    if (strategyResponse.status === 'success' && 
        typeof strategyResponse.data !== 'undefined'){
      setTechSrc(strategyResponse.data.data)
      setTechCol(strategyResponse.data.attributes)
    }
  }, [strategyResponse])

  useEffect(() => {
    if(Init){
      var myChart = echarts.getInstanceByDom(chartRef.current)
      if(techCol.length > 0){
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
    
        // wordcloud data
        // an example
        // var dataList = Object.entries(dataSrc)
        // console.log(dataList)
        // var wc = new Js2WordCloud(wordCloudRef.current);
        // const wc = new Js2WordCloud();
        // static global:
    
        // first plot the proportion of technologies
        myChart.setOption({
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
    
          ],// series end
        })
      }
    } else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      setInit(true)
    }
  }, [src, strategyResponse])

  return (
    <>
      <ToggleAccordion header={
        <h3 id="Statistics">{"Statistics"}
        </h3>}>
        {<div>
          After filtering by keywords, users can take an overview of statistics and word-clouds about the filtered results.
          The left chart is the proportion of sequencing technologies. The middle chart is the proportion of species.
          <br />
          For the construction of the wordclouds (the right chart), we used the term frequency–inverse document frequency strategy (TF-IDF) to evaluate the word weight in dataset summary and overall design based on TfidfVectorizer.
        </div>}
      </ToggleAccordion>
      <Row>
        <Col>
          <div
            ref={chartRef}
            className="chart"
            //the target DOM container needs height and width
            style={{ height: height, width: width, margin: margin }}></div>
        </Col>
        <Col span={8}>
          <div
            ref={wordCloudRef}
            className='wordCloud'
            style={{ height: '22rem', width: '25rem', margin: '2rem' }}
          >
          </div>
        </Col>
      </Row>
    </>
  )
}

DatasetStatic.defaultProps = {
  height: '30rem',
  width: '55rem',
  margin: '0rem',
}

DatasetStatic.propTypes = {
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default DatasetStatic
