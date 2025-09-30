import type {ReactNode} from 'react';
import Link from '@docusaurus/Link';
import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

type FeatureItem = {
  title: string;
  Svg: React.ComponentType<React.ComponentProps<'svg'>>;
  href: string;
  description: ReactNode;
};

const FeatureList: FeatureItem[] = [
  {
    title: 'Retrieval-Augmented Prediction (RAP)',
    Svg: require('@site/static/workflows/rap.svg').default,
    href: '/guide/docs/workflows/rap/overview/',
    description: (
      <>
        RAP is a workflow that uses a retrieval-augmented agent to answer questions about chemical toxicity.
      </>
    ),
  },
  {
    title: 'Systematic Generation of potential Metabolites (SyGMA)',
    Svg: require('@site/static/workflows/sygma.svg').default,
    href: '/guide/docs/workflows/sygma/overview/',
    description: (
      <>
        SyGMA is a workflow that uses a systematic generation of potential metabolites to answer questions about chemical metabolism.
      </>
    ),
  },
  {
    title: 'Ranking',
    Svg: require('@site/static/workflows/ranking.svg').default,
    href: '/guide/docs/workflows/ranking/overview/',
    description: (
      <>
        Ranking is a workflow that uses a pairwise comparison of chemical endpoints to rank chemicals.
      </>
    ),
  },
];

function Feature({title, Svg, href, description}: FeatureItem) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Link to={href} title={title}>
          <Svg className={styles.featureSvg} role="img" />
        </Link>
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures(): ReactNode {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
