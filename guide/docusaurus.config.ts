import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

// This runs in Node.js - Don't use client-side code here (browser APIs, JSX...)

const config: Config = {
  title: 'ToxIndex User Guide',
  tagline: 'Learn how to use ToxIndex effectively',
  favicon: 'logo/logo-no-background.svg',

  // Future flags, see https://docusaurus.io/docs/api/docusaurus-config#future
  future: {
    v4: true, // Improve compatibility with the upcoming Docusaurus v4
  },

  // Set the production url of your site here
  url: 'https://www.toxindex.com',
  // Set the /<baseUrl>/ pathname under which your site is served
  // For GitHub pages deployment, it is often '/<projectName>/'
  baseUrl: '/guide/',
  // Ensure directory-style URLs (/guide/, /docs/page/) resolve to index.html on static hosting
  trailingSlash: true,

  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  organizationName: 'toxindex',
  projectName: 'toxindex-guide',

  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',

  // Even if you don't use internationalization, you can use this field to set
  // useful metadata like html lang. For example, if your site is Chinese, you
  // may want to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  // Plugins
  plugins: [
    [
      require.resolve('@easyops-cn/docusaurus-search-local'),
      {
        hashed: true,
        indexDocs: true,
        indexBlog: false,
        indexPages: true,
      },
    ],
    [
      '@docusaurus/plugin-content-docs',
      {
        id: 'community',
        path: 'community',
        routeBasePath: 'community',
        sidebarPath: require.resolve('./sidebarsCommunity.js'),
        editUrl: undefined,
        showLastUpdateTime: false,
        showLastUpdateAuthor: false,
      },
    ],
    [
      '@docusaurus/plugin-content-docs',
      {
        id: 'developer',
        path: 'developer',
        routeBasePath: 'developer',
        sidebarPath: require.resolve('./sidebarsDeveloper.js'),
        editUrl: undefined,
        showLastUpdateTime: false,
        showLastUpdateAuthor: false,
      },
    ],
  ],

  presets: [
    [
      'classic',
      {
        docs: {
          path: 'docs',
          routeBasePath: 'docs',
          sidebarPath: './sidebars.ts',
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          // editUrl removed
        },
        blog: {
          showReadingTime: false,
          feedOptions: {
            type: ['rss', 'atom'],
            xslt: true,
          },
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          // editUrl removed
          // Useful options to enforce blogging best practices
          onInlineTags: 'warn',
          onInlineAuthors: 'warn',
          onUntruncatedBlogPosts: 'warn',
        },
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    // Replace with your project's social card
    image: 'workflows/sygma.svg',
    navbar: {
      title: 'ToxIndex',
      logo: {
        alt: 'ToxIndex',
        src: 'logo/logo-no-background.svg',
      },
      items: [
        { type: 'docSidebar', sidebarId: 'tutorialSidebar', position: 'left', label: 'Docs' },
        // { type: 'docSidebar', sidebarId: 'communitySidebar', position: 'left', label: 'Community' },
        { to: '/community/intro', label: 'Community', position: 'left' },
        { to: '/developer/intro', label: 'Developer', position: 'left' },
        { href: 'https://github.com/insilica', label: 'GitHub', position: 'right' },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Resources',
          items: [
            { label: 'Getting Started', to: '/docs/intro' },
          ],
        },
        {
          title: 'Workflows',
          items: [
            { label: 'RAP', href: '/docs/workflows/rap/overview/' },
            { label: 'SyGMA', href: '/docs/workflows/sygma/overview/' },
            { label: 'Ranking', href: '/docs/workflows/ranking/overview/' },
          ],
        },
        // {
        //   title: 'Company',
        //   items: [
        //     { label: 'GitHub', href: 'https://github.com/toxindex/toxindex' },
        //   ],
        // },
      ],
      copyright: `© ${new Date().getFullYear()} ToxIndex by Insilica`,
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
