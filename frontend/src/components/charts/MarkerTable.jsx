import React, { useRef, useEffect, useState, useImperativeHandle } from 'react'
import PropTypes from 'prop-types'
import {
  DownloadOutlined,
  EyeOutlined,
  FileTextOutlined,
  LinkOutlined,
  SearchOutlined,
} from '@ant-design/icons'
import strokeColor from '../theme/strokeColor'
import H5adLoader from '../utils/H5adLoader'
import Highlighter from 'react-highlight-words'
import { Button, Input, Space, Table, Tooltip, Row, Col, Select, Spin, Progress } from 'antd'
import * as echarts from 'echarts'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
  TooltipComponent,
  VisualMapComponent,
} from 'echarts/components'
import { LineChart, LinesChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import Loading from './Loading'
import loadingTips from './LoadingTip'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
  VisualMapComponent,
  LinesChart,
])

const { enterLoading, quitLoading } = Loading

const MarkerTable = ({ file, setCompLoad, onRef, title, height, width, margin, progress }) => {
  const commandRef = useRef('')  // current command
  const [action, setAction] = useState(0)
  const [Init, setInit] = useState(false)  // is Inited
  const [loadings, setLoadings] = useState([false, false]) // loading for [marker, heatmap]

  const chartRef = useRef(null) // get current DOM container

  const [lfc, setLFC] = useState([])  // logfoldchanges
  const [mn, setMN] = useState([])  // marker names
  const [pa, setPA] = useState([]) // p-value_adj
  const [me, setME] = useState([]) // mean expression
  const [fr, setFR] = useState([]) // fracments
  const [den, setDen] = useState({}) // dendrogram
  const geneNum = 2000
  const [data, setData] = useState([])
  const [arrLen, setArrLen] = useState(0) // Length of 1D lfc, mn and pval
  const [cluLen, setCluLen] = useState(0) // Length of clusters, arrLen / genes
  const [showTop, setShowTop] = useState(10) // Top markers of each cluster to show
  const [currTip, setCurrTip] = useState(loadingTips[0])
  // filter params
  const [searchText, setSearchText] = useState('')
  const [searchedColumn, setSearchedColumn] = useState('')
  const searchInput = useRef(null)
  const [key, setKey] = useState(true)

  // // 滚动条
  // var xData = key? cells : genes
  // let standardLen = 20
  // var isXTooLong = xData.length > standardLen;
  const [xData, setXData] = useState([])
  const standardLen = 20
  const [Genes, setGenes] = useState([])
  const [Cells, setCells] = useState([])

  const handleSearch = (selectedKeys, confirm, dataIndex) => {
    confirm()
    setSearchText(selectedKeys[0])
    setSearchedColumn(dataIndex)
  }
  const handleReset = (clearFilters) => {
    clearFilters()
    setSearchText('')
  }
  const getColumnSearchProps = (dataIndex) => ({
    filterDropdown: ({
      setSelectedKeys,
      selectedKeys,
      confirm,
      clearFilters,
      close,
    }) => (
      <div
        style={{
          padding: 8,
        }}
        onKeyDown={(e) => e.stopPropagation()}>
        <Input
          ref={searchInput}
          placeholder={`Search ${dataIndex}`}
          value={selectedKeys[0]}
          onChange={(e) =>
            setSelectedKeys(e.target.value ? [e.target.value] : [])
          }
          onPressEnter={() => handleSearch(selectedKeys, confirm, dataIndex)}
          style={{
            marginBottom: 8,
            display: 'block',
          }}
        />
        <Space>
          <Button
            type="primary"
            onClick={() => handleSearch(selectedKeys, confirm, dataIndex)}
            icon={<SearchOutlined />}
            size="small"
            style={{
              width: 90,
            }}>
            Search
          </Button>
          <Button
            onClick={() => clearFilters && handleReset(clearFilters)}
            size="small"
            style={{
              width: 90,
            }}>
            Reset
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              confirm({
                closeDropdown: false,
              })
              setSearchText(selectedKeys[0])
              setSearchedColumn(dataIndex)
            }}>
            Filter
          </Button>
          <Button
            type="link"
            size="small"
            onClick={() => {
              close()
            }}>
            close
          </Button>
        </Space>
      </div>
    ),
    filterIcon: (filtered) => (
      <SearchOutlined
        style={{
          color: filtered ? '#1677ff' : undefined,
        }}
      />
    ),
    onFilter: (value, record) =>
      record[dataIndex].toString().toLowerCase().includes(value.toLowerCase()),
    onFilterDropdownOpenChange: (visible) => {
      if (visible) {
        setTimeout(() => searchInput.current?.select(), 100)
      }
    },
    render: (text) =>
      searchedColumn === dataIndex ? (
        <Highlighter
          highlightStyle={{
            backgroundColor: '#ffc069',
            padding: 0,
          }}
          searchWords={[searchText]}
          autoEscape
          textToHighlight={text ? text.toString() : ''}
        />
      ) : (
        text
      ),
  })
  const initColumns = [
    {
      title: 'Gene',
      dataIndex: 'name',
      key: 'name',
      width: '15%',
      fixed: 'left',
      ...getColumnSearchProps('name'),
    },
    {
      title: 'LogFoldChange',
      dataIndex: 'lfc',
      key: 'lfc',
      sorter: (a, b) => a.key - b.key,
      sortDirections: ['descend', 'ascend'],
    },
    {
      title: 'P-value_adj',
      dataIndex: 'pval',
      key: 'pval',
      sorter: (a, b) => a.key - b.key,
      sortDirections: ['descend', 'ascend'],
    },
    {
      title: 'Cluster',
      dataIndex: 'clu',
      key: 'clu',
      filters: [],
      width: '30%',
      filterSearch: true,
      onFilter: (value, record) => record['clu'] == value,
      sorter: (a, b) => a.key - b.key,
      sortDirections: ['descend', 'ascend'],
    },
  ]
  const [dyColumns, setDyColumns] = useState(initColumns)

  const ScH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['marker', 'den'])
          setMN(h5info.names)
          setLFC(h5info.lfc)
          setPA(h5info.pval)
          setME(h5info.expr)
          setFR(h5info.frac)
          setDen({
            'icoord': h5info.icoord,
            'dcoord': h5info.dcoord,
            'ivl': h5info.ivl,
            'leaves': h5info.leaves,
            'color_list': h5info.color_list,
            'leaves_color_list': h5info.leaves_color_list
          })
          setArrLen(h5info.names.length)
          setCluLen(h5info.names.length / geneNum)
          resolve({
            names: h5info.names,
            lfc: h5info.lfc,
            pval: h5info.pval,
          })
        } catch (error) {
          reject(error)
        }
      }
      reader.onloadend = () => {
        //console.log("Load markers finished!.")
      }
      reader.onerror = (error) => {
        reject(error)
      }
      reader.readAsArrayBuffer(file)
    })

  }

  const toggleAnno = (command) => {
    commandRef.current = command
    setAction(action + 1)
  }

  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }))

  useEffect(() => {
    if (Init) {
      var myChart = echarts.getInstanceByDom(chartRef.current)
      /*...deal with data */
      if (commandRef.current === 'Reload') {
        enterLoading(0, setLoadings)
        enterLoading(1, setLoadings)
        setCurrTip(loadingTips[1])
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['MarkerTable'] = true
          return newCompLoad
        })
        ScH5adLoader(file).then(() => {
          console.log("All data reloaded in MarkerTable.")
          toggleAnno("Upload")
        }).catch(error => {
          console.error('Error fetching blob in MarkerTable:', error)
        })
      }
      if (commandRef.current === "Upload") {
        enterLoading(0, setLoadings)
        enterLoading(1, setLoadings)
        setCurrTip(loadingTips[2])
        let _data = []
        for (let i = 0; i < cluLen; i++) {
          let _arr = mn.slice(i * geneNum, i * geneNum + showTop)  // Top marker for each cluster
          let _obj = _arr.map((item, index) => ({
            'name': item,
            'lfc': lfc[i * geneNum + index],
            'pval': pa[i * geneNum + index],
            'expr': me[i * geneNum + index],
            'frac': fr[i * geneNum + index],
            'clu': i.toString(),
          }))
          _data = [..._data, ..._obj]
        }
        let _clu = [...new Set(_data.map((item) => item['clu']))]
        let _columns = [
          {
            title: 'Gene',
            dataIndex: 'name',
            key: 'name',
            width: '15%',
            fixed: 'left',
            ...getColumnSearchProps('name'),
          },
          {
            title: 'LogFoldChange',
            dataIndex: 'lfc',
            key: 'lfc',
            sorter: (a, b) => a['lfc'] - b['lfc'],
            sortDirections: ['descend', 'ascend'],
          },
          {
            title: 'P-value_adj',
            dataIndex: 'pval',
            key: 'pval',
            sorter: (a, b) => a['pval'] - b['pval'],
            sortDirections: ['descend', 'ascend'],
          },
          {
            title: 'Cluster',
            dataIndex: 'clu',
            key: 'clu',
            filters: _clu.map((item) => ({ 'text': item, 'value': item })),
            filterMode: 'tree',
            filterSearch: true,
            width: '30%',
            onFilter: (value, record) => record['clu'] === value,
            sorter: (a, b) => a['clu'] - b['clu'],
            sortDirections: ['descend', 'ascend'],
          }
        ]
        setData(_data)
        setDyColumns(_columns)
        quitLoading(0, setLoadings)

        console.log(_data)

        // set dendrogram
        let den_data = []
        for (let i = 0; i < den.color_list?.length; i++) {
          den_data.push({
            coords: den.icoord.slice(4 * i, 4 * i + 4).map((item, index) => {
              return [item, den.dcoord[4 * i + index]]
            }),
          })
        }
        // set heatmap in dendrogram order
        let cells = den.ivl
        let genes_unsorted = _data.map((item) => {
          return item.name
        })
        let genes = []

        for (let i = 0; i < cells.length; i++) {
          genes = [
            ...genes,
            ...genes_unsorted.slice(
              den.leaves[i] * showTop,
              (den.leaves[i] + 1) * showTop
            ),
          ]
        }
        genes = genes.reverse()
        _data = _data.sort((a, b) => cells.indexOf(a.clu) - cells.indexOf(b.clu))
        _data = _data.sort((a, b) => genes.indexOf(a.name) - genes.indexOf(b.name))
        console.log(_data)
        console.log(genes)
        let source_test = _data.map((item) => {
          return [item.clu, item.name, item.expr, item.frac]
        })
        console.log(source_test)
        console.log(cells)

        setGenes(genes)
        setCells(cells)

        if (key)
          setXData(cells)
        else
          setXData(genes)

        var isXTooLong = xData.length > standardLen

        console.log(key)
        console.log(isXTooLong)
        console.log(xData)
        myChart.setOption({
          grid: [
            {
              top: '12%',
              left: 25,
              right: '15%',
              bottom: '10%',
              axisLine: {
                lineStyle: {
                  color: '#fff',
                },
              },
              axisPointer: {
                lineStyle: {
                  color: '#ffbd67',
                },
              },
              containLabel: true,
            },
            {
              top: '4%',
              left: 88,
              right: '17%',
              height: '8%',
              show: false,
            },
          ],
          yAxis: [
            {
              gridIndex: 0,
              name: 'Top Markers',
              nameLocation: 'middle',
              nameGap: 70,
              nameRotate: 90,
              nameTextStyle: {
                fontWeight: 'bold',
                fontSize: 16,
              },
              type: 'category',
              splitArea: {
                show: true,
                interval: showTop,
              },
              axisLabel: {
                formatter: function (value, index) {
                  // margin: x label
                  // total: 60 labels
                  // var margin = Math.round(showTop/10)+1
                  // return index % (margin*2) === 0 ? value : '';
                  let margin = Math.round(genes.length / 70)
                  // console.log(margin)
                  return index % margin == 0 ? value : ''
                },
                rotate: 0,
              },
              axisTick: {
                show: true,
                interval: function (index) {
                  // 让刻度线间隔与 axisLabel 一致
                  let margin = Math.round(genes.length / 70)
                  return index % margin === 0
                },
              },
            },
            {
              gridIndex: 1,
              min: 0,
              max: Math.max.apply(null, den.dcoord),
            }
          ],
          xAxis: [
            {
              gridIndex: 0,
              data: cells,
              axisLabel: {
                rotate: 0,
              }
            }, {
              gridIndex: 1,
              min: 0,
              max: Math.max.apply(null, den.icoord),
            }
          ],
          dataset: {
            dimensions: ['cell', 'gene', 'expr', 'frac'],
            source: _data.map((item) => {
              return [item.clu, item.name, item.expr, item.frac]
            }),
          },
          series: [
            {
              type: 'scatter',
              xAxisIndex: 0,
              yAxisIndex: 0,
              encode: {
                x: 'cell',
                y: 'gene',
                tooltip: [0, 1, 2, 3],
              },
              symbolSize: (val) => {
                return val[3] * 18
              },
              emphasis: {
                itemStyle: {
                  borderColor: '#fff',
                  borderWidth: 1,
                },
              },
            },
            {
              type: 'lines',
              xAxisIndex: 1,
              yAxisIndex: 1,
              coordinateSystem: 'cartesian2d',
              polyline: true,
              data: den_data,
              lineStyle: { color: '#gray', width: 2 },
              drawPoint: true,
              silent: false,
            },],
          dataZoom: isXTooLong ? [
            {
              type: 'slider',
              xAxisIndex: 0,
              start: 0,
              end: (standardLen / xData.length) * 100,
              show: true,
            }
          ] : [
            {
              type: 'slider',
              show: false,
              xAxisIndex: 0,
              start: 0,
              end: 100,
            },
          ],
        })
        console.log(myChart.dataZoom)
        quitLoading(1, setLoadings)
      }
    }
    else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      enterLoading(0, setLoadings)
      enterLoading(1, setLoadings)
      setCurrTip(loadingTips[0])
      myChart.setOption({
        textStyle: {
          fontFamily: 'Arial',
        },
        tooltip: {},
        title: [
          {
            text: title,
            left: '40%',
            textStyle: {
              fontSize: 24,
            },
          },
        ],
        grid: [
          {
            top: '12%',
            left: '20',
            right: '15%',
            bottom: '5%',
            // containLabel: false,
            axisLine: {
              lineStyle: {
                color: '#fff',
              },
            },
            axisPointer: {
              lineStyle: {
                color: '#ffbd67',
              },
            },
          },
          {
            top: '4%',
            left: '88',
            right: '21%',
            height: '8%',
            show: false,
          },
        ],
        yAxis: [
          {
            gridIndex: 0,
            name: 'Top Markers',
            nameLocation: 'middle',
            nameGap: 70,
            nameRotate: 90,
            nameTextStyle: {
              fontWeight: 'bold',
              fontSize: 16,
              padding: 1,
            },
            type: 'category',
            splitArea: {
              show: true,
              interval: 4,
            },
            axisTick: {
              alignWithLabel: true,
              // show: true
            },
            axisLabel: {
              interval: 0,
              rotate: 0,
              margin: 10,
              fontSize: 10,
              // show: true,
              // overflow: 'break', // 标签可以超出
            },
          },
          {
            gridIndex: 1,
            type: 'value',
            show: false,
            axisTick: {
              show: false,
            },
            axisLabel: {
              show: false,
            },
            axisLine: {
              show: false,
            },
          },
        ],
        xAxis: [
          {
            gridIndex: 0,
            name: 'Cell Types',
            nameLocation: 'middle',
            nameGap: 30,
            nameRotate: 0,
            nameTextStyle: {
              fontWeight: 'bold',
              fontSize: 16,
              padding: 1,
            },
            type: 'category',
            splitArea: {
              show: true,
            },
            axisTick: {
              alignWithLabel: true,
              // show: true,
            },
            axisLabel: {
              interval: 0,
              fontSize: 14,
              rotate: 0,
              margin: 12,
              // show: true,
              // overflow: 'break', // 标签可以超出
            },
          },
          {
            gridIndex: 1,
            type: 'value',
            show: false,
            axisTick: {
              show: false,
            },
            axisLabel: {
              show: false,
            },
            axisLine: {
              show: false,
            },
          },
        ],
        visualMap: [
          {
            left: 'right',
            top: '23%',
            orient: 'vertical',
            calculable: true,
            dimension: 'frac',
            itemWidth: 20,
            itemHeight: 120,
            inRange: {
              // 设置映射的颜色范围为蓝色到红色之间
              symbolSize: [3, 18],
              color: ["#104e8b", "#ffdab9", "#8b0a50"],
            },
            text: ['Fraction of cells\n in group(%)'],
            textGap: 20,
            textStyle: {
              fontSize: 14,
            },
            padding: 8,
            precision: 1,
            label: {
              formatter: '{value} %',
            },
            min: 0,
            max: 1,
            seriesIndex: 0,
            handleIcon: 'circle',
            handleSize: '100%',
          },
          {
            left: 'right',
            bottom: '23%',
            orient: 'vertical',
            calculable: true,
            precision: 1,
            min: 0,
            max: 5,
            itemHeight: 120,
            dimension: 'expr',
            seriesIndex: 0,
            text: ['Mean Expression\nin group'],
            textGap: 13,
            textStyle: {
              fontSize: 14,
            },
            inRange: {
              color: ["#104e8b", "#ffdab9", "#8b0a50"],
            },
            indicatorIcon: 'roundRect',
            indicatorSize: '60%',
          },
        ],
      })
      setInit(true)
    }
  }, [action, title])// 当 dataPath 改变时重新执行 useEffect

  return (
    <>
      <Row>
        <Col span={24} offset={1}>
          <div style={{ marginBottom: 12 }}>
            Show{' '}
            <Select
              value={showTop}
              options={[5, 10, 15, 30, 50, 100].map(item => ({ value: item }))}
              onChange={(value) => {
                setShowTop(value)
                setKey(true)
                setXData(Cells)
                toggleAnno("Upload")
              }} />
            {' '} Markers of each Cluster
          </div>
        </Col>
        <Col span={22} offset={1}>
          <Spin
            spinning={loadings[0]}
            size="large"
            tip={
              <div>{currTip}
                <Progress percent={progress} strokeColor={strokeColor} size={[300, 15]} />
              </div>
            }>
            <Table
              columns={dyColumns}
              dataSource={data}
              size="small"
              scroll={{
                y: 1000,
              }}
              pagination={{
                showSizeChanger: true,
                showQuickJumper: true,
                showTotal: (total) => `Total ${total} items`,
                pageSizeOptions: [5, 10, 15, 30, 50, 100],
                position: ['bottomRight'],
              }}
            />
          </Spin>
        </Col>
      </Row>
      <Row>
        <Spin
          spinning={loadings[0]}
          size="large"
          tip={
            <div>{currTip}
              <Progress percent={progress} strokeColor={strokeColor} size={[300, 15]} />
            </div>
          }>
          <Button
            onClick={() => {
              // chart update

              var mychart = echarts.getInstanceByDom(chartRef.current)
              const currentOption = mychart.getOption()

              // dataste : [x, y, value1, value2] -> [y, x, value1, value2]
              console.log(currentOption)
              var source = currentOption.dataset[0].source
              var dimensions = currentOption.dataset[0].dimensions

              // exchange dimensions
              dimensions = [dimensions[1], dimensions[0], dimensions[3], dimensions[2]]

              // exchange source
              for (let i = 0; i < source.length; i++) {
                let temp_data = source[i]
                source[i] = [temp_data[1], temp_data[0], temp_data[3], temp_data[2]]
              }

              console.log(source)
              console.log(dimensions)

              const newXAxis = currentOption.yAxis.map(yAxisItem => ({
                ...yAxisItem,
                data: currentOption.yAxis[0].data,
              }))

              const newYAxis = currentOption.xAxis.map(xAxisItem => ({
                ...xAxisItem,
                data: currentOption.xAxis[0].data,
              }))

              if (key) {
                newXAxis[0].axisLabel = {
                  rotate: 90,
                }
                newXAxis[0].nameRotate = 0
                newYAxis[0].nameRotate = 90
              }
              else {
                newXAxis[0].axisLabel = {
                  rotate: 0,
                }
                newYAxis[0].axisLabel = {
                  rotate: 0,
                }
                newXAxis[0].nameRotate = 0
                newYAxis[0].nameRotate = 90
              }

              // 修改 series 的 encode，确保 x 和 y 轴的映射也被交换
              const newSeries = currentOption.series.map(series => ({
                ...series,
                encode: {
                  x: key ? 'gene' : 'cell',
                  y: key ? 'cell' : 'gene',
                  tooltip: [0, 1, 2, 3],
                },
              }))
              console.log(newSeries)

              if (key) {
                setXData(Genes)
              }
              else {
                setXData(Cells)
              }
              var newDataZoom
              if (key)
                newDataZoom = {
                  type: 'slider',
                  xAxisIndex: 0,
                  start: 0,
                  end: (standardLen / Genes.length) * 100,
                  show: true,
                }
              else
                newDataZoom = []


              console.log(xData)

              const newOption = {
                ...currentOption,
                // xAxis: currentOption.yAxis,
                // yAxis: currentOption.xAxis,
                xAxis: newXAxis,
                yAxis: newYAxis,
                grid: [
                  {
                    top: '12%',
                    left: 20,
                    right: '15%',
                    bottom: '10%',
                    axisLine: {
                      lineStyle: {
                        color: '#fff',
                      },
                    },
                    axisPointer: {
                      lineStyle: {
                        color: '#ffbd67',
                      },
                    },
                    containLabel: true,
                  },
                  {
                    top: '4%',
                    left: 88,
                    right: '17%',
                    height: '8%',
                    show: false,
                  },
                ],
                dataset: {
                  dimensions: dimensions,
                  source: source
                },
                series: newSeries,
                dataZoom: newDataZoom,
              }


              mychart.setOption(newOption, true) // 更新图表

              if (key) {
                setKey(false)
              }
              else {
                setKey(true)
              }

            }}>rotation</Button>
          <div
            ref={chartRef}
            className="chart"
            //the target DOM container needs height and width
            style={{
              height: height,
              width: width,
              margin: margin
            }}></div>
        </Spin>
      </Row>
    </>
  )
}

MarkerTable.defaultProps = {
  title: 'Marker Tables',
  height: '60rem',
  width: '60rem',
  margin: '2rem',
}

MarkerTable.propTypes = {
  file: PropTypes.object,
  setCompLoad: PropTypes.func,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
  progress: PropTypes.number,
}

export default MarkerTable
