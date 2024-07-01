import { useState, useEffect } from 'react'
import contactImg from "../assets/img/contact.png"
import layerImg from "../assets/img/LayerView.png"
import pairImg from "../assets/img/pair.png"
import lassoImg from "../assets/img/lasso.png"
import netImg from "../assets/img/network.png"
import heatImg from "../assets/img/heatmap.png"
import coverImg from "../assets/img/cover-figure.png"
import axios from 'axios'
import { Button, ConfigProvider, Layout, Space } from 'antd'
import { ArrowRightCircle } from 'react-bootstrap-icons'
import 'animate.css'
import TrackVisibility from 'react-on-screen'
import { Card, Row, Col } from 'antd'
import { DeliveredProcedureOutlined, SyncOutlined } from '@ant-design/icons'
import { useLocation, useNavigate } from 'react-router-dom'
import Loading from './charts/Loading'
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
  const contentStyle = {
    height: '160px',
    color: '#fff',
    lineHeight: '50px',
    textAlign: 'center',
    background: '#364d79',
  }


  useEffect(() => {
    let ticker = setInterval(() => {
      tick()
    }, delta)

    return () => {
      clearInterval(ticker)
    }
  }, [text])

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
                        <p>
                          Spatial Transcriptomics analysis requires paired single-cell references for cell-type annotation, cell-proportion inference, and cell co-localization detection.
                          However, lacking of paired single-cell and spatial transcriptomics data hinders the in-depth deciphering of tissue architectures and functions.
                          Existing databases mainly provide preliminary pre-defined annotations, instead of sufficiently precise annotations for cell sub-types or cell states of user’s interests, which may not conform to user expectations for integrative analysis.
                         </p>
                         <p>
                          Pairpot curated currently available paired single-cell and spatially resolved transcriptomics data including nearly 830,000 spots from 80 paired studies, with precise annotations for cell sub-types or cell states.
                          Pairpot performed in-depth analysis including batch effects correction, spatial clustering, markers detection, cell-proportion inference, and cell-cell interaction for the paired data.
                          Moreover, Pairpot proposed three online lasso-based analysis modules and integrated them with pre-analysis results into an easy-to-use webpage. In particular, Pairpot designed a LassoView module to refine the customized domains lassoed by users in millisecond-level response time, a PairView module to infer cell proportions of spots using the user-selected cells from single-cell data in real-time, and a LayerView module to display slices in a 3D hierarchical layout.
                        </p>
                        <Button style={{ margin: 0, fontWeight: 500 }}
                          size='large'
                          type='primary'
                          icon={loading[0] ? <SyncOutlined spin /> : <DeliveredProcedureOutlined />}
                          loading={loading[0]}
                          onClick={() => {
                            enterLoading(0, setLoading)
                            setLoadText("Loading Example...")
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
                              quitLoading(0, setLoading)
                              setLoadText("Example Prepared.")
                              navigate('/browse', { state: state })
                            })
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
                      <img src={coverImg} alt="Header Img" />
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
                style={{ width: '46%', margin: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Layer View</h3>}
                cover={<img src={layerImg} ></img>}>
                <Card.Meta classNames='header' description={"The users can explore spatial clusters, cell-type proportions, and gene expression of multiple slices from a study in 3D hierarchical layouts."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '46%', margin: '2%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Pair View</h3>}
                cover={<img src={pairImg} ></img>}>
                <Card.Meta classNames='header' description={"The users can infer cell-type proportions of spots (right) online using the user-lassoed cells from single-cell data (left) in real time."}></Card.Meta>
              </Card>
            </Row>

            <Row>
              <Card hoverable
                style={{ width: '40%', margin: '2%', marginLeft: '2%', marginRight: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Lasso View</h3>}
                cover={<img src={lassoImg} ></img>}>
                <Card.Meta classNames='header' description={"The users can select the customized domains/cells using lasso tools. These domains/cells can be precisely refined in millisecond-level response time."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '26%', margin: '2%', marginRight: '1%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Network</h3>}
                cover={<img src={netImg} ></img>}>
                <Card.Meta classNames='header' description={"The users can explore the cell-cell communication network both in single-cell and spaital data."}></Card.Meta>
              </Card>
              <Card hoverable
                style={{ width: '25%', margin: '2%', marginRight: '2%', marginLeft: '1%', paddingLeft: '1rem', paddingRight: '1rem' }}
                title={<h3 style={{ marginTop: "1rem", padding: 0 }}> Heatmap</h3>}
                cover={<img src={heatImg} ></img>}>
                <Card.Meta classNames='header' description={"The users can explore L-R pairs, single-cell marker genes and spatial variable genes via heatmap in versatile format."}></Card.Meta>
              </Card>
            </Row>


          </section>
        </Content>

        <Footer>
          <Row>

          </Row>
          <Row>
            <Col span={4}>
              <Space direction='horizontal'>
                <h3>Contact us
                  {' '}
                  <img src={contactImg} style={{ width: '20%', height: '20%' }}></img>
                </h3>

              </Space>
            </Col>
            <Col offset={1}>
              <p>
                *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                Medicine, Nankai University, rrrzhan@nankai.edu.cn
              </p>
              <p>
                *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@nankai.edu.cn
              </p>
            </Col>
          </Row>
        </Footer>
      </Layout>
    </ConfigProvider>
  )
}


export default Help
