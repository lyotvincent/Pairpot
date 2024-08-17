import React, { useState, useEffect, View } from 'react'
// import TrackVisibility from 'react-on-screen'
import 'bootstrap/dist/css/bootstrap.min.css'
// import { ReactComponent as Female } from '../assets/img/female.svg'
// import { ReactComponent as FemaleOrgans } from '../assets/img/female-organs.svg'
//import SideMenu from './utils/SideMenu'
import Species from './metas/Species'
import Technologies from './metas/Technologies'
import Tissues from './metas/Tissues'
import Topics from './metas/Topics'
import axios from 'axios'
import { useQuery } from 'react-query'
import { Card, Layout, theme, Tag, Col, Divider, ConfigProvider } from 'antd'
import logoFig from "../assets/img/mylogo.png"
import TagCollapse from './utils/TagCollapse'
// import loadingTips from './charts/LoadingTip'
import Loading from './charts/Loading'
import {
  DatabaseOutlined,
  ProfileOutlined,
  OneToOneOutlined,
} from '@ant-design/icons'
import DatasetTab from './utils/DatasetTab'
import Search from './utils/Search'
import SampleTab from './utils/SampleTab'
const { Content, Footer } = Layout
const { enterLoading } = Loading

export const Database = () => {
  //const [collapsed, setCollapsed] = useState(false)
  const [srcMeta, setSrcMeta] = useState({}) // get the meta data from childrens
  const [src, setSrc] = useState([1, 2, 3])
  const [selectedTags, setSelectedTags] = React.useState(null);

  const handleChange = (tag, checked, label) => {
    if(checked){
      setSelectedTags(tag)
      Rerender(label, tag)
    } else{
      setSelectedTags(null)
      Rerender('all', tag)
    }
  };
  //const [selectedKey, setSelectedKey] = useState('datasets')
  //const [visible, setVisible] = useState(true)
  // const onClickTopMenu = (e) => {
  //   setSelectedKey(e.key)
  // }

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

  const filterSpecies = (species) => {
    return (response.data.data.filter((item) => item[3]?.split(";").includes(species)))
  }

  const filterTissues = (tissues) => {
    return (response.data.data.filter((item) => item[4]?.split(";").includes(tissues)))
  }

  const filterTechs = (species) => {
    return (response.data.data.filter((item) => item[12]?.split(";").includes(species)))
  }

  const filterDiseases = (disease) => {
    return (response.data.data.filter((item) => item[15]?.toLowerCase().includes(disease.toLowerCase())))
  }

  useEffect(() => {
    if (response.status === 'success' && typeof response.data !== 'undefined') {
      setSrc(response)
    }
  }, [response.data])

  const topMenuItems = [
    {
      label: 'Datasets',
      key: 'datasets',
      icon: <DatabaseOutlined />,
    },
    {
      label: 'Samples',
      key: 'samples',
      icon: <ProfileOutlined />,
    },
    {
      label: 'Pairs',
      key: 'pairs',
      icon: <OneToOneOutlined />,
    },
  ]
  const { token } = theme.useToken()
  const TabRef = React.createRef()

  const Rerender = (label, item) => {
    // console.log(label)
    // console.log(item)
    // console.log(src)
    let newSrcData
    if (label === "tissues") {
      newSrcData = filterTissues(item)
      TabRef.current?.GraphConfig(`{tissues: \"${item}\"}`)
    }
    else if (label === "species") {
      newSrcData = filterSpecies(item)
      TabRef.current?.GraphConfig(`{species: \"${item}\"}`)
    }
    else if (label === "technologies") {
      newSrcData = filterTechs(item)
      TabRef.current?.GraphConfig(`{technologies: \"${item}\"}`)
    }
    else if (label === "disease") {
      newSrcData = filterDiseases(item)
    } else if (label === 'all') {
      newSrcData = response.data.data
    } else{
      newSrcData = src.data.data
    }
    setSrc(prevSrc => ({
      ...prevSrc, data: {
        attributes: prevSrc.data.attributes,
        data: newSrcData,
        // add lable and item
        label: label,
        item: item,
      }
    }))
    TabRef.current?.Fresh(true)
    enterLoading(0, TabRef.current?.Loading)
  }

  const species = []

  return (
    <Layout
      style={{
        minHeight: '100vh',
      }}>
      {/* <Sider
        width={200}
        breakpoint="lg"
        collapsible
        collapsed={collapsed}
        onCollapse={(value) => setCollapsed(value)}>
        <SideMenu items={srcMeta} collapsed={collapsed} />
      </Sider> */}
      <Layout>
        <Content
          style={{
            padding: '1rem 1rem',
          }}>
          <Col span={6} offset={9} id='Search'>
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
          <Col span={16} offset={4}> <Search /> </Col>
          <Col span={16} offset={4}>
          <Divider style={{marginTop:16, marginBottom: 14}}/>
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Species.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'species')
                }}>{tag}</Tag.CheckableTag>))}
              threshold={5}
              prefix={<b style={{ margin: 3 }}>Species:</b>}
            />
            {/* <Flex gap={6} wrap justify='left'>
            Species:
            <Tag color="magenta"
              onClick={() => {
                let newSrcData = filterSpecies("Homo sapiens")
                setSrc(prevSrc => ({
                  ...prevSrc, data: {
                    attributes: prevSrc.data.attributes,
                    data: newSrcData
                  }
                }))
                TabRef.current?.Fresh(true)
                enterLoading(0, TabRef.current?.Loading)
              }}>
              Homo sapiens
            </Tag>
            <Tag color="red" onClick={() => {
              Rerender("species", 'Mus musculus')
            }}>Mus musculus</Tag>
            <Tag color="volcano" onClick={() => {
              Rerender("species", 'Oryctolagus cuniculus')
            }}>
              Oryctolagus cuniculus
            </Tag>
            <Tag color="orange" onClick={() => {
              Rerender("species", 'Danio rerio')
            }}>
              Danio rerio
            </Tag>
            <Tag color="gold" onClick={() => {
              Rerender("species", 'Gallus gallus')
            }}>
              Gallus gallus
            </Tag>
          </Flex> */}
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Tissues.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'tissues')
                }}>{tag}</Tag.CheckableTag>))}
              threshold={10}
              prefix={<b style={{ margin: 3 }}>Tissues:</b>}
            />
            {/* <Flex gap={6} wrap justify='left'>
            Tissues:
            <Tag color="magenta"
              onClick={() => {
                Rerender('tissues', 'Brain')
              }}>
              Brain
            </Tag>
            <Tag color="red"
              onClick={() => {
                Rerender('tissues',"Skin")
              }}>
              Skin
            </Tag>
            <Tag color="volcano"
              onClick={() => {
                Rerender('tissues',"Breast")
              }}>
              Breast
            </Tag>
            <Tag color="orange"
              onClick={() => {
                Rerender('tissues',"Kidney")
              }}>
              Kidney
            </Tag>
            <Tag color="gold"
              onClick={() => {
                Rerender('tissues',"Heart")
              }}>
              Heart
            </Tag>
            <Tag color="lime"
              onClick={() => {
                Rerender('tissues',"Lung")
              }}>
              Lung
            </Tag>
            <Tag color="green"
              onClick={() => {
                Rerender('tissues',"Muscle")
              }}>
              Muscle
            </Tag>
            <Tag color="cyan"
              onClick={() => {
                Rerender('tissues',"Prostate")
              }}>
              Prostate
            </Tag>
            <Tag color="blue"
              onClick={() => {
                Rerender('tissues',"Overian")
              }}>
              Overian
            </Tag>
            <Tag color="geekblue"
              onClick={() => {
                Rerender('tissues',"Pancreas")
              }}>
              Pancreas
            </Tag>
            <Tag color="purple"
              onClick={() => {
                Rerender('tissues',"Spinal Cord")
              }}>
              Spinal Cord
            </Tag>
          </Flex> */}
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Technologies.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'technologies')
                }}>{tag}</Tag.CheckableTag>))}
              threshold={6}
              prefix={<b style={{ margin: 3 }}>Platforms:</b>}
            />
            {/* <Flex gap={6} wrap justify='left'>
              Platforms:
              <Tag color="magenta"
                onClick={() => {
                  Rerender('technologies', "10x Visium")
                }}>
                10x Visium
              </Tag>
              <Tag color="red"
                onClick={() => {
                  Rerender('technologies', "scRNA")
                }}>
                scRNA-seq
              </Tag>
              <Tag color="volcano"
                onClick={() => {
                  let newSrcData = filterTechs("snRNA")
                  setSrc(prevSrc => ({
                    ...prevSrc, data: {
                      attributes: prevSrc.data.attributes,
                      data: newSrcData
                    }
                  }))
                  TabRef.current?.Fresh(true)
                  enterLoading(0, TabRef.current?.Loading)
                }}>
                snRNA-seq
              </Tag>
              <Tag color="orange"
                onClick={() => {
                  Rerender('technologies', "Stereo-Seq")
                }}>
                Stereo-Seq
              </Tag>
              <Tag color="gold"
                onClick={() => {
                  Rerender('technologies', "MERFISH")
                }}>
                MERFISH
              </Tag>
              <Tag color="lime"
                onClick={() => {
                  Rerender('technologies', "Spatial Transcriptomics")
                }}>
                Spatial Transcriptomics
              </Tag>
            </Flex> */}
          </Col>
          <Col span={18} offset={3}>
          <TagCollapse
              tags={Topics.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'disease')
                }}>{tag}</Tag.CheckableTag >))}
              threshold={6}
              prefix={<b style={{ margin: 3 }}>Topics:</b>}
            />

            {/* <Flex gap={6} wrap justify='left'>
              Topics:
              <Tag color="magenta"
                onClick={() => {
                  Rerender("disease", "cancer")
                }}>
                Cancer
              </Tag>
              <Tag color="red"
                onClick={() => {
                  Rerender("disease", "injury")
                }}>
                Injury
              </Tag>
              <Tag color="volcano"
                onClick={() => {
                  Rerender("disease", "Alzheimer")
                }}>
                Alzheimer
              </Tag>
              <Tag color="orange"
                onClick={() => {
                  Rerender("disease", "melanoma")
                }}>
                Melanoma
              </Tag>
              <Tag color="gold"
                onClick={() => {
                  Rerender("disease", "carcinoma")
                }}>
                Carcinoma
              </Tag>
            </Flex> */}
          </Col>
          <Divider style={{marginTop:14, marginBottom: 18}}/>
          <Card>
            <DatasetTab sendData={setSrcMeta} response={src} onRef={TabRef} />
          </Card>
          {/* <div>{JSON.stringify(Array.from(new Set(src.data?.data.map((item) => item[4]))))}</div> */}
        </Content>
        <Footer
          style={{
            textAlign: 'center',
          }}>
          PairPot Copyright © 2024. Centre for Bioinformatics and Intelligent Medicine, Nankai University.
        </Footer>
      </Layout>
    </Layout>
  )
}
