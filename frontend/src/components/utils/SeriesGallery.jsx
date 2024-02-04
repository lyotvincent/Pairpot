import Carousel from 'react-multi-carousel'
import meter1 from '../../assets/img/Homo-sapiens.svg'
import meter2 from '../../assets/img/Mus-musculus.svg'
import meter3 from '../../assets/img/Danio-rerio.svg'
import meter4 from '../../assets/img/Macaca-fascicularis.svg'
import 'react-multi-carousel/lib/styles.css'

export const SeriesGallery = () => {
  const responsive = {
    superLargeDesktop: {
      // the naming can be any, depends on you.
      breakpoint: { max: 4000, min: 3000 },
      items: 5,
    },
    desktop: {
      breakpoint: { max: 3000, min: 1024 },
      items: 3,
    },
    tablet: {
      breakpoint: { max: 1024, min: 464 },
      items: 2,
    },
    mobile: {
      breakpoint: { max: 464, min: 0 },
      items: 1,
    },
  }
  return (
    <Carousel
      responsive={responsive}
      infinite={true}
      className="owl-carousel owl-theme species-slider">
      <div className="item">
        <h5>Homo sapiens</h5>
        <h6> (Human)</h6>
        <img src={meter1} alt="Image" />
      </div>
      <div className="item">
        <h5>Mus musculus</h5>
        <h6> (House mouse)</h6>
        <img src={meter2} alt="Image" />
      </div>
      <div className="item">
        <h5>Danio rerio</h5>
        <h6> (Zebrafish)</h6>
        <img src={meter3} alt="Image" />
      </div>
      <div className="item">
        <h5>Macaca fascicularis</h5>
        <h6> (Cynomolgus monkeys)</h6>
        <img src={meter4} alt="Image" />
      </div>
    </Carousel>
  )
}

export default SeriesGallery
