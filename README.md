# Maxim Papusha's Blog

A blog about computational chemistry and cheminformatics built with [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/).

> **⚠️ IMPORTANT**: To deploy this website, you need to enable GitHub Pages in your repository settings.  
> See [DEPLOYMENT.md](DEPLOYMENT.md) for detailed setup instructions.

## 🚀 Features

- **Material Design**: Modern, responsive theme with light/dark mode
- **Jupyter Notebook Support**: Embed and render Jupyter notebooks directly
- **Marimo Notebooks**: Examples and documentation for reactive Python notebooks
- **Math Support**: LaTeX rendering with MathJax
- **Code Highlighting**: Syntax highlighting for multiple languages
- **Search**: Full-text search functionality
- **Automated Deployment**: GitHub Actions workflow for continuous deployment

## 📚 Content

This blog covers topics including:

- Computational chemistry methods and applications
- Cheminformatics tools and techniques
- Molecular modeling and simulation
- Data analysis workflows
- Interactive notebooks and tutorials

## 🛠️ Local Development

### Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/maxim-papusha/maxim-papusha.github.io.git
   cd maxim-papusha.github.io
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Serve the site locally:
   ```bash
   mkdocs serve
   ```

4. Open your browser to `http://127.0.0.1:8000`

### Build

To build the static site:

```bash
mkdocs build
```

The built site will be in the `site/` directory.

## 📝 Adding Content

### Blog Posts

Add new markdown files to the `docs/blog/` directory and update the navigation in `mkdocs.yml`.

### Jupyter Notebooks

1. Place your `.ipynb` file in `docs/notebooks/`
2. Add it to the navigation in `mkdocs.yml`
3. The notebook will be automatically rendered

### Marimo Notebooks

Marimo notebooks are reactive Python notebooks. To include them:

1. Export as HTML: `marimo export html notebook.py > notebook.html`
2. Or document the Python code in a markdown file
3. See `docs/notebooks/marimo-example.md` for an example

## 🔄 Deployment

The site is automatically deployed to GitHub Pages when changes are pushed to the `main` branch using GitHub Actions.

**⚠️ SETUP REQUIRED**: Before the site can be deployed, you must enable GitHub Pages in your repository settings. See [DEPLOYMENT.md](DEPLOYMENT.md) for complete setup instructions.

### Quick Setup

1. Go to **Settings → Pages** in your repository
2. Set **Source** to **"GitHub Actions"**
3. Push a commit or re-run the workflow
4. Your site will be available at `https://maxim-papusha.github.io/`

For detailed troubleshooting and technical information, see [DEPLOYMENT.md](DEPLOYMENT.md).

## 📄 License

This blog content is personal work by Maxim Papusha.

## 🔗 Links

- [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
- [MkDocs](https://www.mkdocs.org/)
- [mkdocs-jupyter](https://github.com/danielfrg/mkdocs-jupyter)
- [Marimo](https://marimo.io/)
