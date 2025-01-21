import { useState, useEffect } from 'react'
import contactImg from "../assets/img/contact.png"
import layerImg from "../assets/img/browse/layer-view.png"
import pairImg from "../assets/img/browse/pair-view.svg"
import lassoImg from "../assets/img/browse/lasso-view.png"
import netImg from "../assets/img/browse/network.svg"
import heatImg from "../assets/img/browse/heatmap.svg"
import coverImg from "../assets/img/cover-figure.png"
import PairpotFooter from './Footer'
//import { Circle } from 'react-bootstrap-icons'
import axios from 'axios'
import { Button, ConfigProvider, Divider, Layout, Space, Statistic } from 'antd'
import { ArrowRightCircle, Cpu } from 'react-bootstrap-icons'
import 'animate.css'
import TrackVisibility from 'react-on-screen'
import { Card, Row, Col } from 'antd'
import { DeliveredProcedureOutlined, SyncOutlined } from '@ant-design/icons'
import { useLocation, useNavigate } from 'react-router-dom'
import Loading from './charts/Loading'
import { useMutation, useQuery, QueryCache } from "react-query"
const { enterLoading, quitLoading } = Loading
const { Meta } = Card

const { Content, Footer } = Layout

const Help = () => {
  const [loopNum, setLoopNum] = useState(0)
  const [isDeleting, setIsDeleting] = useState(false)
  const [text, setText] = useState('')
  const [delta, setDelta] = useState(300 - Math.random() * 100)
  const [index, setIndex] = useState(1)
  const toRotate = ['Web Developer', 'Web Designer', 'UI/UX Designer']
  const period = 2000
  const navigate = useNavigate()
  const [loading, setLoading] = useState([])
  const [loadText, setLoadText] = useState("Try an example >>>")
  const [allData, setAllData] = useState([])
  const [att, setAtt] = useState([])
  const [countAttr, setCountAttr] = useState([])
  const [countData, setCountData] = useState([])
  const contentStyle = {
    height: '160px',
    color: '#fff',
    lineHeight: '50px',
    textAlign: 'center',
    background: '#364d79',
  }

  const response = useQuery({
    queryKey: ['db'],
    queryFn: () => axios.get('/api/datasets').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })

  useEffect(() => {
    // let ticker = setInterval(() => {
    //   tick()
    // }, delta)

    if (response.status == 'success' && response.data != 'undefined') {

      // console.log(response.data)
      setAllData(response.data.data)
      setAtt(response.data.attributes)

      // spots num
      let spots_index = att.indexOf("spots")
      let spots_num = allData.reduce((accumulator, current) => {
        return accumulator + (current[spots_index] || 0)
      }, 0)

      // paired num
      let paired_num = allData.filter(item => {
        return item[1].startsWith("STDS") && item[26] !== "-1" && item[26] !== "0"
      }).length

      // others
      let species_index = att.indexOf("species")
      const uniqueSpecies = new Set()
      allData.forEach(item => {
        const values = item[species_index]?.split(';') || [] // 拆分字符串，确保不为 undefined
        values
          .map(value => value.trim()) // 去除空格
          .filter(value => value) // 过滤掉空字符串
          .forEach(value => uniqueSpecies.add(value)) // 添加到 Set 中
      })
      const species_num = uniqueSpecies.size

      let tissues_index = att.indexOf("tissues")
      const tissues_num = new Set(
        allData
          .map(item => item[tissues_index])
          .filter(tissue => tissue && tissue.trim()) // 过滤掉空值
      ).size
      // console.log(tissues_num)

      let disease_index = att.indexOf("disease")
      const uniqueDisease = new Set()
      allData.forEach(item => {
        const values = item[disease_index]?.split('|') || [] // 拆分字符串，确保不为 undefined
        values
          .map(value => value.trim()) // 去除空格
          .filter(value => value) // 过滤掉空字符串
          .forEach(value => uniqueDisease.add(value)) // 添加到 Set 中
      })
      const disease_num = uniqueDisease.size

      let technologies_index = att.indexOf("technologies")
      const uniqueTechnologies = new Set()
      allData.forEach(item => {
        const values = item[technologies_index]?.split(';') || [] // 拆分字符串，确保不为 undefined
        values
          .map(value => value.trim()) // 去除空格
          .filter(value => value) // 过滤掉空字符串
          .forEach(value => uniqueTechnologies.add(value)) // 添加到 Set 中
      })
      const technologies_num = uniqueTechnologies.size
      // console.log(uniqueTechnologies)

      var countAttr = ["spots", "pairs", "species", "tissues", "diseases", "technologies"]
      var countData = [1425656, 99, 17, tissues_num, disease_num, 25]
      setCountAttr(countAttr)
      setCountData(countData)
    }

    // return () => {
    //   clearInterval(ticker)
    // }
  }, [response.data])

  const tick = () => {
    let i = loopNum % toRotate.length
    let fullText = toRotate[i]
    let updatedText = isDeleting
      ? fullText.substring(0, text.length - 1)
      : fullText.substring(0, text.length + 1)

    setText(updatedText)

    if (isDeleting) {
      setDelta((prevDelta) => prevDelta / 2)
    }

    if (!isDeleting && updatedText === fullText) {
      setIsDeleting(true)
      setIndex((prevIndex) => prevIndex - 1)
      setDelta(period)
    } else if (isDeleting && updatedText === '') {
      setIsDeleting(false)
      setLoopNum(loopNum + 1)
      setIndex(1)
      setDelta(500)
    } else {
      setIndex((prevIndex) => prevIndex + 1)
    }
  }

  return (
    <ConfigProvider theme={{
      components: {
        Card: {
          colorBgContainer: "#fff",
          headerFontSize: 32,
          headerHeight: 48,
          borderRadiusLG: 12,
          colorTextDescription: "rgba(0, 0, 0, 0.6)",
        }
      }
    }}>
      <Layout>
        <Content>
          <section className="banner" id="home">
            <Row>
              <Col span={14} offset={0}>
                <TrackVisibility>
                  {({ isVisible }) => (
                    <div
                      className={
                        isVisible ? 'animate__animated animate__fadeIn' : ''
                      }>
                      <div className='headline'>
                        <h1>
                          {`Welcome to Pairpot!`}
                          {/* <span
                          className="txt-rotate"
                          dataPeriod="1000"
                          data-rotate='[ "Web Developer", "Web Designer", "UI/UX Designer" ]'>
                          <span className="wrap">{text}</span>
                        </span> */}
                        </h1>
                        <h3>
                          {`A Database with Real-time Lasso-Based Analysis Tailored for Paired Single-cell and Spatial Transcriptomics`}
                        </h3>
                        {/* <h5 style={{color:'green'}}>
                        {countAttr.map((attr, i) => (
                            <span key={i}>
                              {`${attr}: ${countData[i]}`} 
                              {i < countAttr.length - 1 && ' '}
                            </span>
                          ))}
                        </h5> */}
                        <p>
                          Paired single-cell and spatial resolved transcriptomics (SRT) data supplement each other, providing in-depth insights into biological processes and disease mechanisms.
                          Pairpot is a database tailored for paired single-cell and SRT data with real-time heuristic analysis.
                          Pairpot curates 99 high-quality pairs including 1,425,656 spots from 299 datasets, and creates the association networks.
                          Pairpot constructs the curated pairs by integrating multiple slices and establishing potential associations between single-cell and SRT data.
                          <br />
                          <br />
                          On this basis, Pairpot adopts semi-supervised learning that enables real-time heuristic analysis for SC-SP pairs where Lasso-View refines the user-selected SRT domains within milliseconds, Pair-View infers cell proportions of spots based on user-selected cell types in real-time, and Layer-View displays SRT slices using a 3D hierarchical layout.
                          In summary, Pairpot provides valuable database resources and powerful heuristic analysis that empower researchers to streamline single-cell and SRT integrative analysis in real-time.
                        </p>
                        <Button style={{ margin: 0, fontWeight: 500 }}
                          size='large'
                          type='primary'
                          icon={loading[0] ? <SyncOutlined spin /> : <DeliveredProcedureOutlined />}
                          loading={loading[0]}
                          onClick={() => {
                            enterLoading(0, setLoading)
                            setLoadText("Loading Example...")
                            // axios({
                            //   method: 'GET',
                            //   url: '/api/example',
                            //   params: {
                            //     id: "STDS0000235"
                            //   },
                            // }).then((response) => {
                            //   let dataCol = response.data.attributes
                            //   let spitem = response.data.data[0]
                            //   let values = Object.fromEntries(
                            //     dataCol.map((k, i) => [k, spitem[i]])
                            //   )
                            //   let scitem = response.data.data[1]
                            //   let state = {
                            //     st: values
                            //   }
                            //   if (typeof scitem !== 'undefined') {
                            //     let scvalues = Object.fromEntries(
                            //       dataCol.map((k, i) => [k, scitem[i]])
                            //     )
                            //     state['sc'] = scvalues
                            //   }

                            // })
                            quitLoading(0, setLoading)
                            setLoadText("Example Prepared.")
                            navigate('/browse')
                          }
                          }
                        >{loadText}</Button>
                      </div>
                    </div>
                  )}
                </TrackVisibility>
              </Col>
              <Col span={9} offset={1}>
                <TrackVisibility>
                  {({ isVisible }) => (
                    <div
                      className={
                        isVisible ? 'animate__animated animate__zoomIn' : ''
                      }
                      style={{
                        margin: "1rem",
                        marginRight: '2rem',
                      }}
                    >
                      <img src={coverImg} alt="Header Img" style={{ width: '60%', marginLeft: '6rem' }} />
                      <Row style={{ marginBottom: '0rem', marginTop: '1rem' }} justify="space-evenly">
                        <Col span={6}>
                          <Statistic title="Cells" value={1180181} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4} style={{ marginLeft: '20px' }}>
                          <Statistic title="Studies" value={299} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4} >
                          <Statistic title="Pairs" value={99} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4}>
                          <Statistic title="Species" value={17} valueStyle={{ color: 'green' }} />
                        </Col>
                      </Row>
                      <Row justify="space-evenly">
                        <Col span={6} >
                          <Statistic title="Spots" value={1425656} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4} style={{ marginLeft: '20px' }}>
                          <Statistic title="Tissues" value={55} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4}>
                          <Statistic title="Diseases" value={149} valueStyle={{ color: 'green' }} />
                        </Col>
                        <Col span={4}>
                          <Statistic title="Technologies" value={25} valueStyle={{ color: 'green' }} />
                        </Col>
                      </Row>
                    </div>
                  )}
                </TrackVisibility>
              </Col>
            </Row>
            <br />

            <Row gutter={16}>
              {/* 
                      <Row>
            <h2>Web Function</h2>
            <img src={funcImg} style={{width:'4%'}}></img>
          </Row>
            <Col span={8}>
              <Card
                cover={
                  <img
                    alt="example"
                    src={propImg}
                    width="40%"
                  />

                }
                actions={[
                  <SettingOutlined key="setting" />,
                  <EditOutlined key="edit" />,
                  <EllipsisOutlined key="ellipsis" />,
                ]}
              >
                <Meta
                  avatar={<Avatar src="https://api.dicebear.com/7.x/miniavs/svg?seed=8" />}
                  title="Card title"
                  description="This is the description"
                />
              </Card>
            </Col> */}
              {/* <Col span={6}>
              <Card 
              cover={
                <img
                    alt="database"
                    src={propImg}
                    style={{width:'53%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              title="Database" bordered={false}>
                The datasets of single-cell and spatial transcriptomes.
              </Card>
            </Col>
            <Col span={6}>
              <Card 
               cover={
                <img
                    alt="browser"
                    src={spImg}
                    style={{width:'65%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              
              title="Browser" bordered={false}>
                Some cell data correspondence plots.
              </Card>
            </Col> */}
              {/* <Col span={8}>
              <Card 
              cover={
                <img
                    alt="help"
                    src={helpImg}
                    style={{width:'56%',display:'block',margin:'auto'}}
                  // src="https://gw.alipayobjects.com/zos/rmsportal/JiqGstEfoWAOHiTxclqi.png"
                  />
               }
              title="Help" bordered={false}>
                Get in touch with us for question help.
              </Card>
            </Col> */}
            </Row>

            {/* <br/>
          <Row>
          <h2>Some data</h2>
          <img src={displayImg} style={{width:'4%'}}></img>
          </Row>
          <div className='datashow'>
            <Row>
              <div className='datalist'>
                <Col span={25}>
                <Card title="Browser" bordered={false}>
                  <table>
                    <tr>
                      <th>hhhhh</th>
                      <th>Species</th>
                    </tr>
                    <tr>
                      <td>A spatial transcriptome map of developing maize ear</td>
                      <td>Maize</td>
                    </tr>
                    <tr>
                      <td>
                      Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis
                      </td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>SAW: An efficient and accurate data analysis workflow for Stereo-seq spatial transcriptomics</td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex infla ...</td>
                      <td>Homo sapiens</td>
                    </tr>
                    <tr>
                      <td>Spatial transcriptomics: The effect of consecutive slices data integration on accurate cell ...</td>
                      <td>Homo sapiens</td>
                    </tr>
                    <tr>
                      <td>Probing pathways by which rhynchophylline modifies sleep using spatial transcriptomics (add ...</td>
                      <td>Mus musculus</td>
                    </tr>
                    <tr>
                      <td>Probing pathways by which rhynchophylline modifies sleep using spatial transcriptomics</td>
                      <td>Mus musculus</td>
                    </tr>
                  </table>
                </Card>
              </Col>
              </div>
              <div className='dataprop'></div>
            </Row>
          </div> */}

            <br />
            <Row>
              <h2 style={{ paddingLeft: '1.5rem' }}>
                {"What can Pairpot provide?"}
              </h2>
            </Row>
            <Row>
              <Card hoverable
                style={{ width: '17%', margin: '0.8%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Layer-View</h3>}
                cover={<img src={layerImg} style={{ height: '10rem' }}></img>}
                onClick={() => {
                  const myVal = "Layer-View"
                  navigate('/browse', { state: { myVal } })
                }}>
                <Card.Meta classnames='header' description={
                  // "The users can explore spatial clusters, cell-type proportions, and gene expression of multiple slices from a study in 3D hierarchical layouts."
                  "Explore spatial clusters, cell-type proportions, and gene expression of multiple slices in 3D hierarchical layouts."
                }></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '17%', margin: '0.8%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Pair-View</h3>}
                cover={<img src={pairImg} style={{ height: '10rem' }}></img>}
                onClick={() => {
                  const myVal = "Pair-View"
                  navigate('/browse', { state: { myVal } })
                }}>
                <Card.Meta classnames='header' description={
                  // "The users can infer cell-type proportions of spots (right) online using the user-lassoed cells from single-cell data (left) in real time."
                  "Infer cell-type proportions of spots online using the user-lassoed cells from single-cell data."
                }></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '17%', margin: '0.8%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Lasso-View</h3>}
                cover={<img src={lassoImg} style={{ height: '10rem' }}></img>}
                onClick={() => {
                  const myVal = "Lasso-View"
                  navigate('/browse', { state: { myVal } })
                }}>
                <Card.Meta classnames='header' description={
                  // "The users can select the customized domains/cells using lasso tools. These domains/cells can be precisely refined in millisecond-level response time."
                  "Select the customized domains/cells using lasso tools."
                }></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '17%', margin: '0.8%', marginRight: '1%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Network</h3>}
                cover={<img src={netImg} style={{ height: '10rem' }}></img>}
                onClick={() => {
                  const myVal = "Interactions"
                  navigate('/browse', { state: { myVal } })
                }}>
                <Card.Meta classnames='header' description={
                  "Explore the cell-cell communication network both in single-cell and spaital data."
                }></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '17%', margin: '0.8%', marginRight: '2%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Heatmap</h3>}
                cover={<img src={heatImg} style={{ height: '10rem' }}></img>}
                onClick={() => {
                  const myVal = "Spatial Markers"
                  navigate('/browse', { state: { myVal } })
                }}>
                <Card.Meta classnames='header' description=
                  {
                    // "The users can explore L-R pairs, single-cell marker genes and spatial variable genes via heatmap in versatile format."
                    "Explore L-R pairs, single-cell marker genes and spatial variable genes via heatmap in versatile format."
                  }></Card.Meta>
              </Card>
            </Row>

            <br />
            <Row gutter={16}>
              <div className='headline'>
                <h2>
                  {"News"}
                </h2>
                <Divider />
                <p>
                  <span style={{ color: 'gray' }}>2024.11.4&nbsp;&nbsp;</span>
                  Pairpot has been published on Nucleic Acids Research online.
                </p>

                <p>
                  <span style={{ color: 'gray' }}>2025.1.16&nbsp;&nbsp;</span>
                  Pairpot has updated some new sc-sp paired dataset: 
                  STDS0000242-SCDS0000020, STDS0000244-SCDS0000021,  STDS0000245-SCDS0000022,STDS0000247-SCDS0000023,STDS0000248-SCDS0000024,STDS0000264-SCDS0000033,STDS0000265-SCDS0000034.         </p>

              </div>
            </Row>


          </section>
        </Content>

        <PairpotFooter/>
      </Layout>
    </ConfigProvider>
  )
}


export default Help
